#!/usr/bin/env python3
"""
analyse_backemf.py — diagnose the back-EMF response of the electric
                     potential Phi to the injected velocity perturbation.

Physical question
-----------------
The Lorentz body force decomposes as

    F = sigma( -grad(Phi) x B )  +  sigma( (U.B) B )  -  sigma |B|^2 U

The last term is the Hartmann body damping that sets the laminar U_theta.
The FIRST term encodes the *electric* response: any U-perturbation feeds
into the charge-conservation equation

    div(sigma grad Phi) = div(sigma U x B)

so a velocity perturbation U' produces a perturbation Phi' that, in turn,
modifies J' = sigma(-grad Phi' + U' x B) and thus the body force.

If |Phi'| is non-negligible compared to its scale (U' B L), the decay we
observed in analyse_perturb.py is electromagnetically mediated, not purely
hydrodynamic.  If |Phi'| is dwarfed by U' B L, the perturbation is being
killed purely by body-force damping on the perpendicular components.

This script
-----------
1. Reads `<probes-dir>/Phi` and the U from the same probes directory.
2. Defines the laminar baseline Phi_lam[i] = Phi at the first sample of each
   probe (= the laminar value, since Phi was untouched by perturb_ic.py).
3. Computes Phi'(t) = Phi(t) - Phi_lam at every sample.
4. Computes the analytic-scale estimate

       Phi_scale(t) = RMS(U_z, t) * B * L

   with B = 0.5 T (axial field magnitude) and L = 5 mm (perturbation
   integral scale).  This is the back-EMF amplitude that a single eddy of
   the injected scale would produce instantaneously, ignoring the
   div-grad-Phi smoothing.
5. Plots
       (a) Phi'(t) per probe
       (b) RMS(Phi') vs Phi_scale  with the ratio overlaid
6. Writes a short verdict comparing the two:
       ratio < 0.1   -> back-EMF response is suppressed (Phi smoothing
                        dominates: low-Rm + closed domain + insulating BCs
                        flattens the Phi' field).
       0.1 < ratio < 1.0 -> back-EMF participates in the dynamics
       ratio > 1.0   -> back-EMF dominates (rare; would imply Hartmann
                        body force is being driven through the Phi route)

Usage
-----
    python3 scripts/analyse_backemf.py --probes-dir postProcessing_5pct/probes/0.5

Defaults to postProcessing/probes/0.5 if no flag.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    HAS_MPL = True
except ImportError:
    HAS_MPL = False


# Physical constants
B_FIELD = 0.5            # T, applied axial magnetic field magnitude
L_INTEG = 5.0e-3         # m, perturbation integral scale (matches perturb_ic.py)
TAU_A = 7.4e-3           # s, Hartmann body-force damping time


# -----------------------------------------------------------------------------
# Probe-file parsers (vector and scalar)
# -----------------------------------------------------------------------------
def read_probes_vector(path: Path) -> tuple[np.ndarray, np.ndarray, list[tuple]]:
    """OpenFOAM vector probe time series. Returns (t, U[Nt,Np,3], locs)."""
    locs: list[tuple] = []
    t: list[float] = []
    U: list[list[list[float]]] = []
    paren = re.compile(r"\(\s*([^)]+)\s*\)")
    with path.open() as f:
        for line in f:
            ln = line.strip()
            if ln.startswith("# Probe"):
                m = re.search(r"\(\s*([\-\d\.eE+ ]+)\s*\)", ln)
                if m:
                    locs.append(tuple(float(x) for x in m.group(1).split()))
                continue
            if not ln or ln.startswith("#"):
                continue
            head, _, rest = ln.partition(" ")
            try:
                tcurr = float(head)
            except ValueError:
                continue
            tuples = paren.findall(rest)
            row = [
                [float(x) for x in tup.replace(",", " ").split()] for tup in tuples
            ]
            t.append(tcurr)
            U.append(row)
    if not t:
        raise RuntimeError(f"{path}: no rows parsed")
    return np.asarray(t), np.asarray(U), locs


def read_probes_scalar(path: Path) -> tuple[np.ndarray, np.ndarray, list[tuple]]:
    """OpenFOAM scalar probe time series. Returns (t, Phi[Nt,Np], locs)."""
    locs: list[tuple] = []
    t: list[float] = []
    Phi: list[list[float]] = []
    with path.open() as f:
        for line in f:
            ln = line.strip()
            if ln.startswith("# Probe"):
                m = re.search(r"\(\s*([\-\d\.eE+ ]+)\s*\)", ln)
                if m:
                    locs.append(tuple(float(x) for x in m.group(1).split()))
                continue
            if not ln or ln.startswith("#"):
                continue
            parts = ln.split()
            try:
                tcurr = float(parts[0])
            except (ValueError, IndexError):
                continue
            try:
                vals = [float(x) for x in parts[1:]]
            except ValueError:
                continue
            t.append(tcurr)
            Phi.append(vals)
    if not t:
        raise RuntimeError(f"{path}: no rows parsed")
    return np.asarray(t), np.asarray(Phi), locs


# -----------------------------------------------------------------------------
# Driver
# -----------------------------------------------------------------------------
def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--probes-dir",
        default="postProcessing/probes/0.5",
        help="probes start time directory",
    )
    p.add_argument(
        "--outdir",
        default=None,
        help="output directory (default: <probes-dir>/../../analyse_backemf/)",
    )
    p.add_argument("--start", type=float, default=None)
    p.add_argument("--end", type=float, default=None)
    args = p.parse_args()

    pdir = Path(args.probes_dir)
    Phi_path = pdir / "Phi"
    U_path = pdir / "U"
    if not Phi_path.is_file():
        print(f"FATAL: {Phi_path} not found.", file=sys.stderr)
        return 1
    if not U_path.is_file():
        print(f"FATAL: {U_path} not found.", file=sys.stderr)
        return 1

    if args.outdir is None:
        # postProcessing_<label>/probes/0.5  ->  postProcessing_<label>/analyse_backemf
        outdir = pdir.parent.parent / "analyse_backemf"
    else:
        outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Read both files
    print(f"Reading {Phi_path}")
    t_phi, Phi, locs = read_probes_scalar(Phi_path)
    print(f"  parsed {len(t_phi)} rows, {Phi.shape[1]} probes")

    print(f"Reading {U_path}")
    t_u, U, _ = read_probes_vector(U_path)
    print(f"  parsed {len(t_u)} rows, {U.shape[1]} probes")

    # Align time grids -- they should be identical (both written timeStep)
    if len(t_phi) != len(t_u) or not np.allclose(t_phi, t_u):
        n = min(len(t_phi), len(t_u))
        t_phi = t_phi[:n]
        Phi = Phi[:n]
        t_u = t_u[:n]
        U = U[:n]
        print(f"  truncated to common length {n}")
    t = t_phi

    if args.start is not None:
        m = t >= args.start
        t = t[m]; Phi = Phi[m]; U = U[m]
    if args.end is not None:
        m = t <= args.end
        t = t[m]; Phi = Phi[m]; U = U[m]

    tau = t - t[0]
    Nt, Np = Phi.shape

    # Phi'(t) = Phi(t) - Phi(t=t0).  At t0, Phi is laminar; perturbation
    # was injected into U only.
    Phi_lam = Phi[0:1]                   # (1, Np)
    Phi_p = Phi - Phi_lam                # (Nt, Np)

    # Analytic back-EMF scale: Phi_scale = RMS(U_z) * B * L
    Uz = U[..., 2]
    rms_Uz = np.sqrt(np.mean(Uz ** 2, axis=1))                 # (Nt,)
    Phi_scale = rms_Uz * B_FIELD * L_INTEG                     # (Nt,)

    # RMS(Phi') across probes
    rms_Phi_p = np.sqrt(np.mean(Phi_p ** 2, axis=1))           # (Nt,)

    # Ratio measures the strength of the realised back-EMF compared to the
    # naive eddy estimate.
    ratio = np.divide(
        rms_Phi_p, Phi_scale,
        out=np.zeros_like(rms_Phi_p),
        where=Phi_scale > 0,
    )

    # Comparison reference: laminar Phi RMS at the probes
    Phi_lam_rms = np.sqrt(np.mean(Phi_lam ** 2))

    # Mean ratio over the body of the run (skip first 2 settling steps)
    body = slice(2, None)
    ratio_med = float(np.median(ratio[body])) if Nt > 2 else float("nan")
    ratio_max = float(ratio[body].max())     if Nt > 2 else float("nan")

    # Verdict
    if ratio_med < 0.1:
        verdict = "WEAK back-EMF (decay is essentially hydrodynamic)"
    elif ratio_med < 1.0:
        verdict = "MODERATE back-EMF (electromagnetic feedback participates)"
    else:
        verdict = "STRONG back-EMF (Phi route dominates)"

    print()
    print("=== Back-EMF diagnostics ===")
    print(f"  N_probes                      = {Np}")
    print(f"  time window                   = {t[0]:.4f}..{t[-1]:.4f} s  "
          f"(span {tau[-1]*1e3:.2f} ms = {tau[-1]/TAU_A:.2f} tau_A)")
    print(f"  laminar Phi RMS (at probes)   = {Phi_lam_rms:.4e} V")
    print(f"  initial RMS Phi' (induced)    = {rms_Phi_p[0]:.4e} V")
    print(f"  peak    RMS Phi' (induced)    = {rms_Phi_p.max():.4e} V "
          f"at t-t0 = {tau[rms_Phi_p.argmax()]*1e3:.2f} ms")
    print(f"  final   RMS Phi'              = {rms_Phi_p[-1]:.4e} V")
    print(f"  RMS Phi' / Phi_lam (peak)     = "
          f"{rms_Phi_p.max()/Phi_lam_rms*100:.3f}%")
    print()
    print(f"  back-EMF scale  U'_z B L      = {Phi_scale[0]:.4e} V (initial)")
    print(f"  realised / scale (median)     = {ratio_med:.3f}")
    print(f"  realised / scale (max)        = {ratio_max:.3f}")
    print()
    print(f"  VERDICT: {verdict}")

    # ------------------------------------------------------------------
    # Output text
    # ------------------------------------------------------------------
    vfile = outdir / "verdict_backemf.txt"
    with vfile.open("w") as f:
        f.write(
            f"Back-EMF diagnostic, debug case\n"
            f"  injection time              : t = {t[0]:.4f} s\n"
            f"  observation window          : "
            f"{tau[-1]*1e3:.2f} ms ({tau[-1]/TAU_A:.2f} tau_A)\n"
            f"  laminar Phi RMS             : {Phi_lam_rms:.4e} V\n"
            f"  initial induced RMS Phi'    : {rms_Phi_p[0]:.4e} V\n"
            f"  peak induced RMS Phi'       : {rms_Phi_p.max():.4e} V "
            f"(at t-t0 = {tau[rms_Phi_p.argmax()]*1e3:.2f} ms)\n"
            f"  final induced RMS Phi'      : {rms_Phi_p[-1]:.4e} V\n"
            f"  peak Phi' / laminar Phi     : "
            f"{rms_Phi_p.max()/Phi_lam_rms*100:.3f}%\n"
            f"  back-EMF scale U'_z B L     : {Phi_scale[0]:.4e} V (init)\n"
            f"  realised/scale (median)     : {ratio_med:.3f}\n"
            f"  realised/scale (max)        : {ratio_max:.3f}\n"
            f"  VERDICT                     : {verdict}\n"
        )
    print(f"\nWrote {vfile}")

    if not HAS_MPL:
        print("matplotlib not installed; skipping plots")
        return 0

    # Plot 1: Phi'(t) per probe
    fig, ax = plt.subplots(figsize=(8, 5))
    for i, loc in enumerate(locs):
        ax.plot(tau * 1e3, Phi_p[:, i] * 1e6, label=f"probe {i}: {loc}")
    ax.axhline(0.0, color="k", lw=0.5)
    ax.set_xlabel("t - t0  [ms]")
    ax.set_ylabel(r"$\Phi'$  [$\mu$V]")
    ax.set_title(r"$\Phi'(t) = \Phi(t) - \Phi_{\rm lam}$ at the five probes")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(outdir / "Phi_prime_per_probe.png", dpi=140)
    plt.close(fig)
    print(f"Wrote {outdir/'Phi_prime_per_probe.png'}")

    # Plot 2: RMS Phi' vs analytic eddy scale, plus ratio
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

    ax1.semilogy(tau * 1e3, rms_Phi_p * 1e6, "o-", label=r"RMS $\Phi'$ "
                 r"(realised)")
    ax1.semilogy(tau * 1e3, Phi_scale * 1e6, "r--",
                 label=r"$U'_z\,B\,L$ (eddy scale)")
    ax1.set_ylabel(r"voltage [$\mu$V]")
    ax1.set_title(f"Back-EMF magnitude: {verdict}")
    ax1.legend(fontsize=10, loc="best")
    ax1.grid(alpha=0.3, which="both")

    ax2.plot(tau * 1e3, ratio, "o-")
    ax2.axhline(0.1, color="g", lw=0.5, label="weak/moderate threshold")
    ax2.axhline(1.0, color="r", lw=0.5, label="moderate/strong threshold")
    ax2.set_xlabel("t - t0  [ms]")
    ax2.set_ylabel("realised RMS / eddy scale")
    ax2.set_ylim(bottom=0.0)
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(outdir / "backEMF_ratio.png", dpi=140)
    plt.close(fig)
    print(f"Wrote {outdir/'backEMF_ratio.png'}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
