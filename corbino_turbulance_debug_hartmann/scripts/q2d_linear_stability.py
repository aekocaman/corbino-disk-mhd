#!/usr/bin/env python3
"""
q2d_linear_stability.py — Q2D linear stability of the laminar Corbino flow.

"""

from __future__ import annotations
from pathlib import Path

import numpy as np
import scipy.linalg

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt



RHO    = 6440.0           # kg / m^3
NU     = 3.73e-7          # m^2 / s
SIGMA  = 3.46e6           # S / m
B0     = 0.5              # T

H      = 0.010            # m   gap height
R_IN   = 0.020            # m
R_OUT  = 0.100            # m

MU_DYN = RHO * NU
HA     = (H / 2.0) * np.sqrt(SIGMA * B0**2 / MU_DYN)
TAU_A  = RHO / (SIGMA * B0**2)
TAU_H  = HA * TAU_A


#  Chebyshev differentiation matrix

def cheb_diff(N: int):
    """N+1 Chebyshev points on [-1, 1] and the differentiation matrix.

    Trefethen, Spectral Methods in MATLAB, p. 58.
    Convention: x[0] = +1 (xi=+1 -> r=r_out), x[N] = -1 (xi=-1 -> r=r_in).
    """
    if N == 0:
        return np.zeros((1, 1)), np.array([1.0])
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.r_[2, np.ones(N - 1), 2] * (-1.0) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1.0 / c) / (dX + np.eye(N + 1))
    D = D - np.diag(D.sum(axis=1))
    return D, x



#  Eigenproblem assembly with null-space-projected BCs

def assemble_eigenproblem(m: int, U_inf: float, N: int = 32):
"""Build (A_reduced, B_reduced) for the m-th azimuthal mode.

    The full (N+1) x (N+1) matrices are projected onto the null space
    of the 4-row constraint matrix C that enforces psi = psi' = 0 at
    both walls.  The reduced problem has dimension (N-3) x (N-3) and
    contains only the physically admissible modes.
    """
    Dxi, xi = cheb_diff(N)
    L = R_OUT - R_IN
    r = 0.5 * (R_IN + R_OUT) + 0.5 * L * xi
    Dr = (2.0 / L) * Dxi
    D2r = Dr @ Dr

    Rinv  = np.diag(1.0 / r)
    Rinv2 = np.diag(1.0 / r ** 2)

    Lm  = D2r + Rinv @ Dr - m**2 * Rinv2
    Lm2 = Lm @ Lm

    K = U_inf * 0.5 * (R_IN + R_OUT)
    UoverR = K / r ** 2
    UoverR_diag = np.diag(UoverR)

    A_full = Lm.astype(complex)
    B_full = (NU * Lm2
              - Lm / TAU_H
              - 1j * m * (UoverR_diag @ Lm)).astype(complex)

    # ---- constraint matrix C: 4 rows for psi = psi' = 0 at walls ----
    C = np.zeros((4, N + 1))
    C[0, 0]   = 1.0          # psi(r_out) = 0   (xi = +1, node 0)
    C[1, :]   = Dr[0, :]     # psi'(r_out) = 0
    C[2, :]   = Dr[N, :]     # psi'(r_in)  = 0
    C[3, N]   = 1.0          # psi(r_in)  = 0

    # ---- null space basis Z of C (shape: (N+1, N-3)) ----
    Z = scipy.linalg.null_space(C)
    # Sanity check: Z should have shape (N+1, N-3)
    assert Z.shape == (N + 1, N - 3), (
        f"Null space shape mismatch: got {Z.shape}, expected "
        f"({N + 1}, {N - 3}).  Constraint matrix C may be ill-conditioned."
    )

    # ---- projection ----
    Zh = Z.conj().T
    A_red = Zh @ A_full @ Z
    B_red = Zh @ B_full @ Z
    return A_red, B_red


def linear_stability_q2d(m: int, U_inf: float, N: int = 32) -> float:
    """Return Re(sigma_max) for mode m at U_inf, in 1/s."""
    A, B = assemble_eigenproblem(m, U_inf, N=N)
    eigvals = scipy.linalg.eigvals(B, A)
    eigvals = eigvals[np.isfinite(eigvals)]
    if len(eigvals) == 0:
        return -np.inf
    return float(np.max(np.real(eigvals)))





def main():
    print("Physical setup")
    print(f"  Ha     = {HA:.3f}")
    print(f"  tau_A  = {TAU_A * 1e3:.4f} ms")
    print(f"  tau_H  = {TAU_H * 1e3:.2f} ms")
    print(f"  -1/Ha  = {-1.0 / HA:+.5f}  "
          f"(analytical Q2D prediction for sigma*tau_A)")
    print()

    outdir = Path("results")
    outdir.mkdir(exist_ok=True)

    modes = list(range(1, 9))
    ratios = [0.5, 0.75, 1.0, 1.19, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0]

    sigma_table = np.zeros((len(modes), len(ratios)))

    print(f"{'Re_a/Ha^2':>10} | {'U_inf [cm/s]':>13} | "
          f"{'sigma_tA range over m':>30}")
    print("-" * 62)
    for j, ratio in enumerate(ratios):
        Re_a = ratio * HA**2
        U_inf = Re_a * NU / (R_OUT - R_IN)
        for i, m in enumerate(modes):
            s_real = linear_stability_q2d(m, U_inf, N=32)
            sigma_table[i, j] = s_real * TAU_A
        col = sigma_table[:, j]
        print(f"{ratio:10.3f} | {U_inf*100:13.3f} | "
              f"[{col.min():+.5f}, {col.max():+.5f}]")
    print()

    target = -1.0 / HA
    deviation = np.abs(sigma_table - target).max()
    rel_dev = deviation / abs(target)
    print(f"Reference -1/Ha = {target:+.5f}.")
    print(f"Maximum |sigma_max*tau_A - (-1/Ha)| across sweep: "
          f"{deviation:.2e}  (relative {rel_dev:.2e})")
    print()

    # --- Marginal stability: find where sigma_max*tau_A crosses 0 ---
    sigma_per_ratio = sigma_table.max(axis=0)
    crit = None
    if (sigma_per_ratio > 0).any():
        idx_first_pos = np.where(sigma_per_ratio > 0)[0][0]
        if idx_first_pos > 0:
            r0 = ratios[idx_first_pos - 1]
            r1 = ratios[idx_first_pos]
            s0 = sigma_per_ratio[idx_first_pos - 1]
            s1 = sigma_per_ratio[idx_first_pos]
            crit = r0 + (r1 - r0) * (-s0) / (s1 - s0)
            print(f"sigma_max*tau_A crosses 0 between Re_a/Ha^2 = "
                  f"{r0:.2f} and {r1:.2f};")
            print(f"linear interpolation gives critical "
                  f"(Re_a/Ha^2)_c ~ {crit:.3f}")
        else:
            print("sigma_max*tau_A is positive at the lowest sweep point; "
                  "expand to lower Re_a/Ha^2.")
    else:
        last_ratio = ratios[-1]
        last_sigma = sigma_per_ratio[-1]
        prev_ratio = ratios[-2]
        prev_sigma = sigma_per_ratio[-2]
        slope = (last_sigma - prev_sigma) / (last_ratio - prev_ratio)
        if slope > 0:
            crit = last_ratio + (-last_sigma) / slope
            print(f"sigma_max*tau_A negative across the whole sweep; "
                  f"linear extrapolation gives sigma_max=0 at "
                  f"Re_a/Ha^2 ~ {crit:.2f}")
        else:
            print("sigma_max*tau_A negative across the whole sweep with "
                  "non-monotonic trend; no marginal stability detected.")
    print()


    # Plot

    fig, axs = plt.subplots(1, 2, figsize=(11, 4.5))

    # -- panel (a): sigma_max vs Re_a/Ha^2 --
    ax = axs[0]
    sigma_max_per_ratio = sigma_table.max(axis=0)

    # Reference lines first (so data plots on top)
    ax.axhline(0.0, color="C3", lw=1.0, ls=":",
               label=r"marginal stability $\sigma\tau_A=0$")
    ax.axhline(target, color="0.5", lw=1.0, ls="--",
               label=f"Hermitian limit $-1/Ha = {target:+.4f}$")
    if crit is not None and ratios[0] <= crit <= ratios[-1]:
        ax.axvline(crit, color="C2", lw=0.8, ls=":")
        ax.scatter([crit], [0.0], marker="x", s=110, color="C2",
                   linewidths=2.0, zorder=5,
                   label=fr"$(Re_a/Ha^2)_c \approx {crit:.2f}$")

    # Numerical data on top
    ax.plot(ratios, sigma_max_per_ratio, "o-",
            color="#0d4f8b", lw=1.6, ms=7, zorder=4,
            label=r"$\sigma_{\max}\tau_A$ (max over $m$)")

    ax.set_xlabel(r"$Re_a / Ha^2$")
    ax.set_ylabel(r"$\sigma_{\max}\,\tau_A$")
    ax.set_title(r"(a)  Q2D growth rate vs operating point")
    margin = 0.005
    ax.set_ylim(min(sigma_max_per_ratio.min(), target) - margin,
                sigma_max_per_ratio.max() + margin)
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(loc="upper left", fontsize=9, framealpha=0.95)


    ax = axs[1]
    extent = [ratios[0] - 0.1, ratios[-1] + 0.1,
              modes[0] - 0.5,  modes[-1] + 0.5]
    im = ax.imshow(sigma_table, aspect="auto", origin="lower",
                   extent=extent, cmap="viridis")
    cb = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
    cb.set_label(r"$\sigma\tau_A$")
    ax.set_xlabel(r"$Re_a / Ha^2$")
    ax.set_ylabel(r"azimuthal mode $m$")
    ax.set_title(r"(b)  $\sigma\tau_A(m, Re_a/Ha^2)$")
    ax.set_yticks(modes)

    fig.tight_layout()
    fig.savefig(outdir / "q2d_linstab.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {outdir / 'q2d_linstab.png'}")


    # Summary markdown

    summary = outdir / "q2d_linstab.md"
    with open(summary, "w") as f:
        f.write("# Q2D Linear Stability of the Laminar Corbino Flow\n\n")
        f.write(f"- Hartmann number:       $Ha = {HA:.2f}$\n")
        f.write(f"- Hartmann time:         $\\tau_H = {TAU_H*1e3:.2f}$ ms\n")
        f.write(f"- Analytical prediction: "
                f"$\\sigma\\tau_A \\approx -1/Ha = {target:+.5f}$\n\n")
        f.write("## Computed `sigma_max * tau_A`\n\n")
        f.write("| $Re_a/Ha^2$ | $U_\\infty$ (cm/s) | "
                "$\\sigma_{\\max}\\tau_A$ (max over $m$) | "
                "$|\\sigma_{\\max}\\tau_A - (-1/Ha)|$ |\n")
        f.write("|------------:|------------------:|"
                "-----------------------------------:|"
                "----------------------------------:|\n")
        for j, ratio in enumerate(ratios):
            Re_a = ratio * HA**2
            U_inf = Re_a * NU / (R_OUT - R_IN)
            sigma_max = sigma_table[:, j].max()
            dev = abs(sigma_max - target)
            f.write(
                f"| {ratio:.2f} | {U_inf*100:.3f} | "
                f"${sigma_max:+.5f}$ | "
                f"${dev:.2e}$ |\n"
            )

    print(f"Wrote {summary}")


if __name__ == "__main__":
    main()
