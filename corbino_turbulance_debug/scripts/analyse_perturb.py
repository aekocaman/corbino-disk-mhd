#!/usr/bin/env python3


from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

# Try matplotlib; fall back to text-only output if it's not installed.
try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    HAS_MPL = True
except ImportError:
    HAS_MPL = False


TAU_A = 7.4e-3  # s, Lorentz/Hartmann damping time at our parameters



# Probe-file parser (OpenFOAM ASCII vector probe format)

# Format:
#   # Probe 0 (x y z)
#   # Probe 1 ...
#   #       Probe         0          1          ...
#   #         Time
#       t1   (Ux Uy Uz) (Ux Uy Uz) ...
#
def read_probes(path: Path) -> tuple[np.ndarray, np.ndarray, list[tuple]]:
    """Return (t, U, locs) where U has shape (Nt, Nprobe, 3)."""
    locs: list[tuple] = []
    t: list[float] = []
    U: list[list[list[float]]] = []
    line_re = re.compile(r"\(\s*([^)]+)\s*\)")
    with path.open() as f:
        for line in f:
            line = line.strip()
            if line.startswith("# Probe"):
                m = re.search(r"\(\s*([\-\d\.eE+ ]+)\s*\)", line)
                if m:
                    locs.append(tuple(float(x) for x in m.group(1).split()))
                continue
            if not line or line.startswith("#"):
                continue
            # Numerical row: time, then one (a b c) per probe
            head, _, rest = line.partition(" ")
            try:
                tcurr = float(head)
            except ValueError:
                continue
            tuples = line_re.findall(rest)
            row = [
                [float(x) for x in tup.replace(",", " ").split()]
                for tup in tuples
            ]
            t.append(tcurr)
            U.append(row)
    if not t:
        raise RuntimeError(f"{path}: no numerical rows parsed")
    t_arr = np.asarray(t)
    U_arr = np.asarray(U)
    if U_arr.ndim != 3 or U_arr.shape[2] != 3:
        raise RuntimeError(
            f"{path}: parsed shape {U_arr.shape}, expected (Nt, Nprobe, 3)"
        )
    return t_arr, U_arr, locs



# Diagnostics

def project_theta(U: np.ndarray, locs: list[tuple]) -> np.ndarray:
    """Project U_xy at each probe onto local theta-hat. Returns (Nt, Nprobe)."""
    out = np.zeros((U.shape[0], U.shape[1]))
    for i, (x, y, _z) in enumerate(locs):
        r = np.hypot(x, y)
        if r < 1e-12:
            continue
        # theta-hat in cartesian: (-y/r, x/r, 0)
        out[:, i] = (-y / r) * U[:, i, 0] + (x / r) * U[:, i, 1]
    return out


def fit_log_linear(t: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    """Fit y = A * exp(sigma * t).  Returns (sigma, A, R^2)."""
    mask = (y > 0) & np.isfinite(y)
    if mask.sum() < 3:
        return float("nan"), float("nan"), 0.0
    tt = t[mask]
    yy = np.log(y[mask])
    # Standard linear regression
    n = len(tt)
    sx = tt.sum()
    sy = yy.sum()
    sxx = (tt * tt).sum()
    sxy = (tt * yy).sum()
    denom = n * sxx - sx * sx
    if abs(denom) < 1e-30:
        return float("nan"), float("nan"), 0.0
    slope = (n * sxy - sx * sy) / denom
    intercept = (sy - slope * sx) / n
    A = float(np.exp(intercept))
    sigma = float(slope)
    yhat = slope * tt + intercept
    ss_res = ((yy - yhat) ** 2).sum()
    ss_tot = ((yy - yy.mean()) ** 2).sum() + 1e-300
    r2 = float(1.0 - ss_res / ss_tot)
    return sigma, A, r2



# Driver

def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--probes-dir",
        default="postProcessing/probes/0.5",
        help="probes start time directory",
    )
    p.add_argument("--start", type=float, default=None,
                   help="trim to t >= start")
    p.add_argument("--end", type=float, default=None,
                   help="trim to t <= end")
    p.add_argument(
        "--outdir",
        default="postProcessing/analyse_perturb",
        help="output directory",
    )
    args = p.parse_args()

    probes_path = Path(args.probes_dir) / "U"
    if not probes_path.is_file():
        print(f"FATAL: {probes_path} not found.", file=sys.stderr)
        return 1

    print(f"Reading {probes_path}")
    t, U, locs = read_probes(probes_path)
    print(f"  parsed {len(t)} rows, {U.shape[1]} probes")

    if args.start is not None:
        m = t >= args.start
        t, U = t[m], U[m]
    if args.end is not None:
        m = t <= args.end
        t, U = t[m], U[m]
    print(f"  using {len(t)} rows in t in [{t[0]:.4f}, {t[-1]:.4f}]")

    # Time-since-injection (perturbation was injected at t = 0.5)
    tau = t - t[0]

    # Components
    Uz = U[..., 2]                                        # (Nt, Nprobe)
    Utheta = project_theta(U, locs)                       # (Nt, Nprobe)

    # RMS over probes -- our scalar perturbation amplitude
    rms_z = np.sqrt(np.mean(Uz ** 2, axis=1))             # (Nt,)
    rms_theta = np.sqrt(np.mean(Utheta ** 2, axis=1))     # (Nt,) -- baseline

    # Drop the first two samples for the fit -- the first PISO sweep
    # cleans up divergence and produces a transient blip in U_z that's
    # not part of the perturbation evolution.
    fit_start_idx = 2
    sigma, A, r2 = fit_log_linear(tau[fit_start_idx:], rms_z[fit_start_idx:])

    # Verdict
    sigma_tauA = sigma * TAU_A
    if not np.isfinite(sigma):
        verdict = "INCONCLUSIVE (fit failed)"
    elif sigma > 0:
        if sigma_tauA > 0.3:
            verdict = "GROWTH"
        else:
            verdict = "MARGINAL GROWTH"
    else:
        if abs(sigma_tauA) > 0.3:
            verdict = "DECAY"
        else:
            verdict = "MARGINAL DECAY"

    print("\n=== Perturbation evolution diagnostics ===")
    print(f"  Probes:                {U.shape[1]}")
    print(f"  Time window:           {t[0]:.4f} to {t[-1]:.4f} s "
          f"(span {tau[-1]*1e3:.2f} ms = {tau[-1]/TAU_A:.2f} tau_A)")
    print(f"  RMS_z initial:         {rms_z[0]:.4e} m/s")
    print(f"  RMS_z final:           {rms_z[-1]:.4e} m/s")
    print(f"  Ratio final/initial:   {rms_z[-1]/rms_z[0]:.3f}")
    print(f"  RMS_theta (laminar):   {rms_theta.mean():.4e} m/s")
    print(f"  Initial perturbation / laminar:  "
          f"{rms_z[0]/rms_theta.mean()*100:.2f}%")
    print(f"  Final   perturbation / laminar:  "
          f"{rms_z[-1]/rms_theta.mean()*100:.2f}%")
    print()
    print(f"  Exponential fit log(RMS_z) = log(A) + sigma * t (skip first 2)")
    print(f"  sigma  =  {sigma:+.4e}  1/s")
    print(f"  sigma * tau_A =  {sigma_tauA:+.4f}  (dimensionless)")
    print(f"  R^2    =  {r2:.4f}")
    print(f"  e-folding |1/sigma| = {1/abs(sigma)*1e3:.2f} ms"
          if np.isfinite(sigma) and sigma != 0 else "")
    print(f"\n  VERDICT:  {verdict}")

    # ------------------------------------------------------------------
    # Outputs
    # ------------------------------------------------------------------
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Verdict text
    verdict_txt = outdir / "verdict.txt"
    with verdict_txt.open("w") as f:
        f.write(
            f"Perturbed-restart diagnostic, debug case\n"
            f"  injection time       : t = {t[0]:.4f} s\n"
            f"  observation window   : {tau[-1]*1e3:.2f} ms "
            f"({tau[-1]/TAU_A:.2f} tau_A)\n"
            f"  initial RMS U_z      : {rms_z[0]:.4e} m/s\n"
            f"  final RMS U_z        : {rms_z[-1]:.4e} m/s\n"
            f"  ratio final/initial  : {rms_z[-1]/rms_z[0]:.3f}\n"
            f"  laminar U_theta RMS  : {rms_theta.mean():.4e} m/s\n"
            f"  perturbation/laminar : "
            f"initial {rms_z[0]/rms_theta.mean()*100:.2f}%, "
            f"final {rms_z[-1]/rms_theta.mean()*100:.2f}%\n"
            f"  exponential rate     : sigma = {sigma:+.4e} 1/s\n"
            f"  sigma * tau_A        : {sigma_tauA:+.4f}\n"
            f"  R^2 of log-linear fit: {r2:.4f}\n"
            f"  e-folding time       : "
            f"{1/abs(sigma)*1e3:.2f} ms\n"
            f"  VERDICT              : {verdict}\n"
        )
    print(f"\nWrote {verdict_txt}")

    if not HAS_MPL:
        print("matplotlib not installed; skipping plots")
        return 0

    # Plot 1: U_z per probe vs t
    fig, ax = plt.subplots(figsize=(8, 5))
    for i, loc in enumerate(locs):
        ax.plot(tau * 1e3, Uz[:, i] * 1e3, label=f"probe {i}: {loc}")
    ax.axhline(0.0, color="k", lw=0.5)
    ax.set_xlabel("t - 0.5 [ms]")
    ax.set_ylabel("$U_z$  [mm/s]")
    ax.set_title("$U_z$ at each probe (zero in laminar; carries perturbation)")
    ax.legend(fontsize=8, loc="best")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(outdir / "Uz_per_probe.png", dpi=140)
    plt.close(fig)
    print(f"Wrote {outdir/'Uz_per_probe.png'}")

    # Plot 2: RMS_z(t) with exponential fit
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.semilogy(tau * 1e3, rms_z * 1e3, "o-", label="RMS $U_z$ (probes)")
    if np.isfinite(sigma):
        tt_fit = tau[fit_start_idx:]
        ax.semilogy(
            tt_fit * 1e3,
            A * np.exp(sigma * tt_fit) * 1e3,
            "r--",
            label=f"fit: sigma = {sigma:+.2e} 1/s ($\\sigma \\tau_A$ = {sigma_tauA:+.3f})",
        )
    ax.set_xlabel("t - 0.5 [ms]")
    ax.set_ylabel("RMS $U_z$ across probes [mm/s]")
    ax.set_title(f"Perturbation evolution: {verdict}")
    ax.legend(fontsize=10)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(outdir / "RMSz_with_fit.png", dpi=140)
    plt.close(fig)
    print(f"Wrote {outdir/'RMSz_with_fit.png'}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
