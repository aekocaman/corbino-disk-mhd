#!/usr/bin/env python3
"""
hartmann_profile.py — Sample U_θ(z) at several radii and compare to
the analytical Hartmann-channel solution.

Analytical solution (steady, parallel plates, uniform body force):

    U(z) = U_∞ · [ 1  −  cosh((z − h/2)/δ_H) / cosh(Ha) ]

where δ_H = sqrt(μ / σ B²) = (h/2) / Ha is the Hartmann boundary-layer
thickness, U_∞ is the local Lorentz/damping balance, and Ha is the
Hartmann number based on h/2.  At each radius U_∞ scales as 1/r
(Corbino current J_r ∝ 1/r), so the family of curves over r tests
both the boundary-layer shape and the radial scaling at once.

For our case:  h = 10 mm,  Ha = 95,  δ_H ≈ 53 μm.

Outputs (saved to the case directory):
  - hartmann_profile_r{NN}cm.csv       
  - hartmann_profile.png               
"""

import sys
import numpy as np
import pyvista as pv
from pathlib import Path


RHO     = 6440.0      # kg/m³
NU      = 3.73e-7     # m²/s
MU      = RHO * NU    # dynamic viscosity, Pa·s
SIGMA   = 3.46e6      # S/m
B0      = 0.5         # T  (axial field magnitude)
H_GAP   = 0.01        # m  (channel height between plates)

# Sampling — sweep over radius
R_LIST  = [0.04, 0.06, 0.08]   # m
PHI     = 0.0                  # rad (azimuth — pick a clean column)
N_PTS   = 400                  # samples across z

# Plot colours — fixed convention
COLOR_SIM = "#d62728"   # red
COLOR_ANA = "#1f77b4"   # blue


def az_unit_vector_at(points: np.ndarray) -> np.ndarray:
    """θ̂ unit vector at each (x,y,z); θ̂ = (-sin φ, cos φ, 0)."""
    x, y = points[:, 0], points[:, 1]
    r_xy = np.hypot(x, y)
    r_xy[r_xy == 0] = 1.0
    return np.column_stack((-y / r_xy, x / r_xy, np.zeros_like(x)))


def hartmann_analytic(z, U_inf, h, Ha):
    """Hartmann profile, normalized to U_inf at the channel centre.

    Uses log-domain ratio  cosh(η)/cosh(Ha) = exp(log_cosh(η) − log_cosh(Ha))
    with  log_cosh(x) = |x| − log 2 + log1p(exp(−2|x|))
    to stay finite for any η.
    """
    delta_H = (0.5 * h) / Ha
    eta = (z - 0.5 * h) / delta_H

    abs_eta = np.abs(eta)
    log_cosh_eta = abs_eta - np.log(2.0) + np.log1p(np.exp(-2.0 * abs_eta))
    log_cosh_Ha  = Ha - np.log(2.0) + np.log1p(np.exp(-2.0 * Ha))
    ratio = np.exp(np.minimum(log_cosh_eta - log_cosh_Ha, 0.0))
    return U_inf * (1.0 - ratio)


def sample_one_radius(internal: pv.UnstructuredGrid, r: float, phi: float):
    """Return (z_line, U_theta_sim) along the vertical at (r, φ)."""
    eps = 1.0e-6
    x0, y0 = r * np.cos(phi), r * np.sin(phi)
    line = internal.sample_over_line(
        (x0, y0, 0.0 + eps),
        (x0, y0, H_GAP - eps),
        resolution=N_PTS - 1,
    )
    if "U" in line.point_data:
        U_line = np.asarray(line.point_data["U"])
    else:
        U_line = np.asarray(line.cell_data["U"])

    pts = np.asarray(line.points)
    z_line = pts[:, 2]
    that = az_unit_vector_at(pts)
    U_theta = np.einsum("ij,ij->i", U_line, that)

    valid = (
        (np.linalg.norm(U_line, axis=1) > 0.0)
        & (z_line >= 0.0) & (z_line <= H_GAP)
    )
    return z_line[valid], U_theta[valid]


def main(case_dir: Path):
    foam_file = case_dir / "case.foam"
    if not foam_file.exists():
        foam_file.touch()

    print(f"[load] {foam_file}")
    reader = pv.OpenFOAMReader(str(foam_file))
    times = reader.time_values
    print(f"[load] using last time = {times[-1]}")
    reader.set_active_time_value(times[-1])
    reader.disable_all_patch_arrays()
    reader.enable_patch_array("internalMesh")
    multi = reader.read()
    internal = multi["internalMesh"]
    print(f"[load] {internal.n_cells} cells, fields = "
          f"{list(internal.cell_data.keys())}")

    Ha = 0.5 * H_GAP * np.sqrt(SIGMA * B0**2 / MU)
    delta_H = (0.5 * H_GAP) / Ha
    print(f"[analytic] Ha = {Ha:.2f}   δ_H = {delta_H*1e6:.2f} μm")
    print()

    # ------------------------------------------------------------------- #
    # Sample, fit, and gather diagnostics for every radius
    # ------------------------------------------------------------------- #
    results = []
    for r in R_LIST:
        z_line, U_theta = sample_one_radius(internal, r, PHI)
        i_mid = int(np.argmin(np.abs(z_line - 0.5 * H_GAP)))
        U_inf = U_theta[i_mid]
        U_th_pred = hartmann_analytic(z_line, U_inf, H_GAP, Ha)

        err = U_theta - U_th_pred
        abs_err_max = float(np.max(np.abs(err)))
        rms_full = float(np.sqrt(np.mean(err**2)))
        bulk = (z_line > 5 * delta_H) & (z_line < H_GAP - 5 * delta_H)
        rms_bulk = float(np.sqrt(np.mean(err[bulk]**2)))

        print(f"--- r = {r*100:.0f} cm "
              f"-----------------------------------------------------")
        print(f"  U_∞                     = {U_inf:+.5f}  m/s")
        print(f"  max |U_θ_sim − U_th|    = {abs_err_max:.3e}  m/s "
              f"({100*abs_err_max/abs(U_inf):.2f} % of U_∞)")
        print(f"  rms (full)              = {rms_full:.3e}  m/s "
              f"({100*rms_full/abs(U_inf):.2f} % of U_∞)")
        print(f"  rms (bulk)              = {rms_bulk:.3e}  m/s "
              f"({100*rms_bulk/abs(U_inf):.2f} % of U_∞)")
        # 1/r expectation: r·U_∞ should be invariant
        print(f"  r · U_∞                 = {r*U_inf:+.4e}  m²/s")
        print()

        results.append({
            "r": r, "z": z_line, "U_sim": U_theta, "U_ana": U_th_pred,
            "U_inf": U_inf, "rms_full": rms_full, "rms_bulk": rms_bulk,
        })

        out_csv = case_dir / f"hartmann_profile_r{int(r*100):02d}cm.csv"
        np.savetxt(
            out_csv,
            np.column_stack([z_line, U_theta, U_th_pred]),
            header="z[m], U_theta_sim[m/s], U_theta_analytic[m/s]",
            delimiter=",", comments="",
        )
        print(f"[save] {out_csv}")

    # consistency: r·U_∞ should be constant
    rU = np.array([res["r"] * res["U_inf"] for res in results])
    print(f"\n[1/r check]  r·U_∞  = {rU}")
    print(f"[1/r check]  spread  = "
          f"{100*(rU.max()-rU.min())/abs(rU.mean()):.2f} %")

    # ------------------------------------------------------------------- #
    # Combined figure: one row per r, two columns (full | near-wall)
    # ------------------------------------------------------------------- #
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        n = len(results)
        fig, axs = plt.subplots(
            n, 2, figsize=(11, 3.2 * n),
            gridspec_kw={"width_ratios": [2, 1]}, sharex="col",
        )
        if n == 1:
            axs = np.array([axs])

        for row, res in enumerate(results):
            ax_full, ax_zoom = axs[row, 0], axs[row, 1]
            r = res["r"]; z = res["z"]
            U_sim = res["U_sim"]; U_ana = res["U_ana"]; U_inf = res["U_inf"]

            ax_full.plot(z * 1e3, U_sim, "o", ms=3,
                         color=COLOR_SIM, label="Simulation")
            ax_full.plot(z * 1e3, U_ana, "-", lw=1.5,
                         color=COLOR_ANA, label="Hartmann Analytic")
            ax_full.axhline(U_inf, color="k", ls=":", alpha=0.4)
            ax_full.set_ylabel(r"$U_\theta$  [m/s]")
            ax_full.set_title(
                f"Hartmann profile, r = {r*100:.0f} cm, Ha = {Ha:.0f}"
            )
            # The profile dips from 0 at the walls to U_∞ < 0 in the bulk;
            # empty space sits above the curve, centred over z = h/2.
            ax_full.legend(loc="upper center", framealpha=0.9)
            ax_full.grid(alpha=0.3)

            wall = z < 5 * delta_H
            ax_zoom.plot(z[wall] * 1e6, U_sim[wall], "o-", ms=4,
                         color=COLOR_SIM, label="Simulation")
            ax_zoom.plot(z[wall] * 1e6, U_ana[wall], "-", lw=1.5,
                         color=COLOR_ANA, label="Hartmann Analytic")
            ax_zoom.axvline(delta_H * 1e6, color="k", ls="--", alpha=0.6,
                            label=fr"$\delta_H$ = {delta_H*1e6:.0f} μm")
            ax_zoom.set_ylabel(r"$U_\theta$  [m/s]")
            ax_zoom.set_title("Near-wall region")
            # Curve descends from (z=0, U_θ=0) to (z≈5δ_H, U_θ≈U_∞);
            # upper-right (high z, U_θ near 0) is empty.
            ax_zoom.legend(fontsize=8, loc="upper right", framealpha=0.9)
            ax_zoom.grid(alpha=0.3)

        axs[-1, 0].set_xlabel("z  [mm]")
        axs[-1, 1].set_xlabel("z  [μm]")

        fig.tight_layout()
        out_png = case_dir / "hartmann_profile.png"
        fig.savefig(out_png, dpi=140)
        print(f"[save] {out_png}")
    except Exception as exc:
        print(f"[plot] skipped ({exc})")


if __name__ == "__main__":
    case = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd()
    main(case.resolve())
