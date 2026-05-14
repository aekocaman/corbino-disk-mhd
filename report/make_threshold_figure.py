#!/usr/bin/env python3
"""
make_threshold_figure.py - Produces the threshold-crossing figure
for consolidated_part_v4.tex (Sec. V.E).

Combines the three empirical data points from the Phase 2 bracket
(Re_a/Ha^2 = 1.19), the Phase B perturbation run (Re_a/Ha^2 = 6.0),
and the Phase C-B perturbation run (Re_a/Ha^2 = 14.6) with the Q2D
linear-stability curve from q2d_linstab.py to make the headline
"factor-of-4" plot.

Run from the report/ directory:
    python3 make_threshold_figure.py

Output:
    threshold_crossing.png   (300 dpi, ready for v4 \\includegraphics)
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

# ---------------------------------------------------------------------
# Empirical data: three operating points, single seed (1729), back-half
# log-linear fit of <U_z^2>^{1/2}(t).  See campaign_log entries.
# ---------------------------------------------------------------------
EMPIRICAL = {
    "Re_a/Ha^2":  np.array([1.19,  6.0,    14.6]),
    "sigma_tA":   np.array([-0.298, -0.167, -0.032]),
    "R2":         np.array([0.81,   0.996,  0.56]),
    "label":      "3D simulation (debug mesh, alpha = 2%, seed 1729)",
    "seed_var":   0.10,            # +/- 1-sigma seed-realisation variance
}

# ---------------------------------------------------------------------
# Q2D linear stability sweep, reproduced from q2d_linstab.py output.
# At each Re_a/Ha^2 the maximum real eigenvalue over m in {1,...,8}.
# Approximate values from the existing q2d_linstab.png.  Replace with
# the exact array from results/q2d_linstab.md if you regenerate.
# ---------------------------------------------------------------------
Q2D_RA_HA2 = np.array(
    [0.5, 0.75, 1.0, 1.19, 1.5, 2.0, 3.0, 3.78, 4.0, 5.0, 7.0, 10.0]
)
Q2D_SIGMA_TA = np.array(
    [-0.0105, -0.0103, -0.0100, -0.0097, -0.0091, -0.0079, -0.0048,
      0.0000,  0.0014,  0.0048,  0.0095,  0.0145]
)

# ---------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------
plt.rcParams.update({
    "font.family":         "sans-serif",
    "axes.linewidth":      0.9,
    "axes.labelsize":     12,
    "axes.titlesize":     12,
    "xtick.labelsize":    10,
    "ytick.labelsize":    10,
    "xtick.direction":    "in",
    "ytick.direction":    "in",
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "legend.fontsize":    10,
    "legend.frameon":     False,
    "mathtext.fontset":   "stix",
})


def main():
    fig, ax = plt.subplots(figsize=(7.0, 4.6))

    # ---------- Q2D linear-stability curve ----------
    ax.plot(Q2D_RA_HA2, Q2D_SIGMA_TA, "-", color="#1f77b4", lw=2.0,
            label="Q2D linear stability (Sommeria-Moreau)", zorder=2)
    ax.fill_between(Q2D_RA_HA2, Q2D_SIGMA_TA, 0,
                    where=(Q2D_SIGMA_TA < 0), alpha=0.07, color="#1f77b4",
                    zorder=1)

    # Q2D threshold marker
    ax.axvline(3.78, color="#1f77b4", ls=":", lw=1.0)
    ax.text(3.78, 0.04,
            r"$(Re_{a}/Ha^{2})_{c,\rm Q2D}\approx 3.78$",
            color="#1f77b4", fontsize=10, ha="center",
            bbox=dict(boxstyle="round,pad=0.2", fc="white",
                      ec="#1f77b4", alpha=0.85))

    # ---------- Empirical 3D points ----------
    yerr = np.full_like(EMPIRICAL["sigma_tA"], EMPIRICAL["seed_var"])
    ax.errorbar(EMPIRICAL["Re_a/Ha^2"], EMPIRICAL["sigma_tA"],
                yerr=yerr, fmt="o", color="#cc3333", ms=10,
                mfc="#cc3333", mec="black", mew=0.9,
                ecolor="#cc3333", elinewidth=1.4, capsize=4,
                label="3D simulation (debug mesh, this work)",
                zorder=4)
    ax.plot(EMPIRICAL["Re_a/Ha^2"], EMPIRICAL["sigma_tA"],
            "-", color="#cc3333", lw=1.4, alpha=0.7, zorder=3)

    # Annotate each point
    for x, y, r2 in zip(EMPIRICAL["Re_a/Ha^2"],
                        EMPIRICAL["sigma_tA"],
                        EMPIRICAL["R2"]):
        offset_y = 0.045 if y > -0.1 else -0.06
        ax.annotate(
            rf"$\sigma\tau_A={y:+.3f}$" + "\n" + rf"($R^{{2}}={r2:.2f}$)",
            (x, y), xytext=(x, y + offset_y),
            fontsize=9, color="#7a1d1d", ha="center",
            arrowprops=dict(arrowstyle="-", color="#cc3333",
                            lw=0.5, alpha=0.6),
        )

    # ---------- Empirical 3D threshold (linear interpolation) ----------
    # Linear interp of last two empirical points to sigma_tA = 0:
    x1, y1 = EMPIRICAL["Re_a/Ha^2"][1], EMPIRICAL["sigma_tA"][1]
    x2, y2 = EMPIRICAL["Re_a/Ha^2"][2], EMPIRICAL["sigma_tA"][2]
    Re3D_c = x1 + (0 - y1) * (x2 - x1) / (y2 - y1)
    ax.axvline(Re3D_c, color="#cc3333", ls="--", lw=1.2)
    ax.text(Re3D_c, 0.085,
            rf"$(Re_{{a}}/Ha^{{2}})_{{c,\rm 3D}}\approx {Re3D_c:.0f}$",
            color="#cc3333", fontsize=10, ha="center",
            bbox=dict(boxstyle="round,pad=0.2", fc="white",
                      ec="#cc3333", alpha=0.95))

    # ---------- Marginal stability line ----------
    ax.axhline(0.0, color="black", lw=0.7, ls="-")
    ax.text(0.55, 0.005, "marginal stability", fontsize=9,
            color="0.4", style="italic")

    # ---------- "Factor of 4" arrow between the two thresholds ----------
    arrow = FancyArrowPatch(
        (3.78, -0.18), (Re3D_c, -0.18),
        arrowstyle="<->", color="0.25", lw=1.4,
        mutation_scale=14,
    )
    ax.add_patch(arrow)
    ax.text((3.78 + Re3D_c) / 2, -0.21,
            rf"factor $\sim {Re3D_c/3.78:.1f}\!\times$"
            r" : 3D damping channels",
            fontsize=10, color="0.25", ha="center", style="italic")

    # ---------- Axes ----------
    ax.set_xscale("log")
    ax.set_xlim(0.4, 22)
    ax.set_ylim(-0.45, 0.10)
    ax.set_xlabel(r"$Re_{a}/Ha^{2}$ (control parameter)")
    ax.set_ylabel(r"dimensionless growth rate $\sigma\,\tau_{A}$")
    ax.set_title("Threshold-crossing campaign: Q2D theory vs 3D simulation")
    ax.grid(True, which="both", alpha=0.25, lw=0.5)

    # PK14 region label (optional context)
    ax.axvspan(0.5, 1.5, alpha=0.08, color="#2ca02c", zorder=0)
    ax.text(0.95, -0.42, "PK14 expt.\nregime",
            fontsize=8.5, color="#2ca02c", ha="center")

    leg = ax.legend(loc="lower right", framealpha=0.95)
    leg.get_frame().set_edgecolor("0.7")

    plt.tight_layout()
    out = "threshold_crossing.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"-> {out}")
    print()
    print("Empirical 3D threshold (linear interpolation of last 2 points):")
    print(f"  (Re_a/Ha^2)_c,3D = {Re3D_c:.2f}")
    print(f"  Q2D threshold    = 3.78")
    print(f"  Ratio            = {Re3D_c/3.78:.2f}")


if __name__ == "__main__":
    main()
