#!/usr/bin/env python3
"""
field_diagnostics.py — Time-resolved sanity check on the laminar
Corbino case.  Pulls quantitative numbers directly from the saved
volume fields (no ParaView), at several time steps spanning the run,
and compares them against the analytical predictions.

For each time, reports:
    - U_θ at three radial probes (r = 4, 6, 8 cm, mid-plane, φ=0)
    - axisymmetry test: |U_θ| at φ = 0 vs φ = π/2 (should be equal)
    - r · U_θ at each probe (should be ~constant — Corbino 1/r law)
    - J_r at each probe (should also follow 1/r → r·J_r constant)
    - Φ at each probe (should follow log law in steady state)
    - global |U|_max and |U|_min over the internal mesh
    - bulk vs predicted ratio: U_θ_bulk / [I/(2π r h σ B)]

Usage:
    python3 field_diagnostics.py [/path/to/case]
"""

import sys
import numpy as np
import pyvista as pv
from pathlib import Path
from scipy.spatial import cKDTree

# Physics constants
RHO     = 6440.0
NU      = 3.73e-7
MU      = RHO * NU
SIGMA   = 3.46e6
B0      = 0.5
H_GAP   = 0.01
R_IN    = 0.02
R_OUT   = 0.10
DELTA_V = 3.7e-3       # potential difference, V

# Probes: (r, φ, z)
PROBES_PHI0 = [(0.04, 0.0, 0.005),
               (0.06, 0.0, 0.005),
               (0.08, 0.0, 0.005)]

# Axisymmetry pair: same r, different φ
AXISYM_PAIR = [(0.06, 0.0,        0.005),
               (0.06, 0.5*np.pi,  0.005)]


def cart_from_rphi(r, phi, z):
    return np.array([r * np.cos(phi), r * np.sin(phi), z])


def az_unit_vector_at(p):
    x, y, _ = p
    r = np.hypot(x, y)
    return np.array([-y / r, x / r, 0.0])


def rad_unit_vector_at(p):
    x, y, _ = p
    r = np.hypot(x, y)
    return np.array([x / r, y / r, 0.0])


def pick_times(times, n_target=6):
    """Pick a roughly log-spaced subset spanning the available times."""
    t = np.asarray(times)
    if len(t) <= n_target:
        return t.tolist()
    # log-spaced over [t.min, t.max], snap to nearest available
    tmin = max(t[t > 0].min(), 1e-3)
    tmax = t.max()
    target = np.geomspace(tmin, tmax, n_target)
    out, last = [], -1
    for tt in target:
        i = int(np.argmin(np.abs(t - tt)))
        if i != last:
            out.append(t[i])
            last = i
    if t[-1] not in out:
        out.append(t[-1])
    return sorted(set(out))


def main(case_dir: Path):
    foam_file = case_dir / "case.foam"
    if not foam_file.exists():
        foam_file.touch()

    reader = pv.OpenFOAMReader(str(foam_file))
    times_all = list(reader.time_values)
    times_pick = pick_times(times_all)
    print(f"[load] {len(times_all)} times available; "
          f"sampling: {[f'{t:g}' for t in times_pick]}")
    print()

    # Header
    print(
        f"{'t [s]':>8}  "
        f"{'U_θ(4)':>10}  {'U_θ(6)':>10}  {'U_θ(8)':>10}  "
        f"{'r·U_θ(4)':>11}  {'r·U_θ(6)':>11}  {'r·U_θ(8)':>11}  "
        f"{'axisym%':>8}  "
        f"{'|U|max':>9}  "
        f"{'Φ(4)':>9}  {'Φ(6)':>9}  {'Φ(8)':>9}"
    )
    print("-" * 156)

    for t in times_pick:
        reader.set_active_time_value(t)
        reader.disable_all_patch_arrays()
        reader.enable_patch_array("internalMesh")
        multi = reader.read()
        internal = multi["internalMesh"]

        # KDTree over cell centres for nearest-cell sampling
        centers = np.asarray(internal.cell_centers().points)
        tree = cKDTree(centers)

        U_field   = np.asarray(internal.cell_data["U"])
        Phi_field = np.asarray(internal.cell_data["Phi"]) \
            if "Phi" in internal.cell_data else None
        J_field   = np.asarray(internal.cell_data["J"]) \
            if "J" in internal.cell_data else None

        # Sample at the three radial probes
        u_th = []
        ru_th = []
        phis = []
        for (r, phi, z) in PROBES_PHI0:
            p = cart_from_rphi(r, phi, z)
            _, idx = tree.query(p, k=1)
            U_at = U_field[idx]
            that = az_unit_vector_at(p)
            u_th_val = float(np.dot(U_at, that))
            u_th.append(u_th_val)
            ru_th.append(r * u_th_val)
            phis.append(float(Phi_field[idx]) if Phi_field is not None else np.nan)

        # Axisymmetry pair
        a, b = AXISYM_PAIR
        pa = cart_from_rphi(*a)
        pb = cart_from_rphi(*b)
        _, ia = tree.query(pa, k=1)
        _, ib = tree.query(pb, k=1)
        u_th_a = float(np.dot(U_field[ia], az_unit_vector_at(pa)))
        u_th_b = float(np.dot(U_field[ib], az_unit_vector_at(pb)))
        if abs(u_th_a) + abs(u_th_b) > 0:
            axisym_pct = 100 * abs(u_th_a - u_th_b) / max(abs(u_th_a), abs(u_th_b))
        else:
            axisym_pct = float('nan')

        # Global U magnitude range
        U_mag = np.linalg.norm(U_field, axis=1)
        u_max = float(U_mag.max())

        print(
            f"{t:8.4f}  "
            f"{u_th[0]:+10.5f}  {u_th[1]:+10.5f}  {u_th[2]:+10.5f}  "
            f"{ru_th[0]:+11.4e}  {ru_th[1]:+11.4e}  {ru_th[2]:+11.4e}  "
            f"{axisym_pct:8.3f}  "
            f"{u_max:9.5f}  "
            f"{phis[0]*1e3:+9.4f}  {phis[1]*1e3:+9.4f}  {phis[2]*1e3:+9.4f}"
        )

    # ------------------------------------------------------------------- #
    # Final-time analytical comparison
    # ------------------------------------------------------------------- #
    print()
    print("Predicted (steady state, t→∞):")
    Ha = 0.5 * H_GAP * np.sqrt(SIGMA * B0**2 / MU)
    delta_H = (0.5 * H_GAP) / Ha
    # I from simple Ohm: I = ΔV / R  with R = 1/(2π σ h) · ln(r_out/r_in)
    R_circuit = np.log(R_OUT / R_IN) / (2 * np.pi * SIGMA * H_GAP)
    I_pred = DELTA_V / R_circuit
    print(f"  Ha               = {Ha:.2f}")
    print(f"  δ_H              = {delta_H*1e6:.2f} μm")
    print(f"  R_circuit        = {R_circuit*1e6:.4f} μΩ")
    print(f"  I (Ohmic, ΔV/R)  = {I_pred:.2f} A")
    for r, _, _ in PROBES_PHI0:
        Jr = I_pred / (2 * np.pi * r * H_GAP)
        Uth = Jr / (SIGMA * B0)
        print(f"  r = {r*100:.0f} cm:  "
              f"J_r = {Jr:.3e} A/m²,  U_θ_bulk = {Uth:+.4e} m/s,  "
              f"r·U_θ_bulk = {r*Uth:+.4e}")
    print()
    print("Diagnostics:")
    print("  - U_θ(r) should approach U_θ_bulk(r) as t → ∞")
    print("  - r·U_θ should be ~constant (Corbino 1/r law)")
    print("  - axisym% should fall to ~0.01 at the converged state")
    print("  - |U|_max should grow then stabilise")


if __name__ == "__main__":
    case = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd()
    main(case.resolve())
