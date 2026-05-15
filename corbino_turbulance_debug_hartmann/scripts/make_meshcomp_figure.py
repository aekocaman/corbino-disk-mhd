#!/usr/bin/env python3
"""
make_meshcomp_figure.py — produce the mesh-sensitivity comparison plot.


"""

from __future__ import annotations
from pathlib import Path
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt



# Data: debug-mesh same-seed bracket (5 amplitudes, seed 1729)

DEBUG = {
    "amp_pct":  np.array([2.0, 5.0, 10.0, 15.0, 30.0]),
    "sigma_tA": np.array([-0.298, -0.310, -0.332, -0.353, -0.404]),
    "R2":       np.array([0.81, 0.84, 0.87, 0.91, 0.95]),
    "tfold":    np.array([24.9, 23.9, 22.3, 21.0, 18.3]),
}




HRES = {
    "amp_pct_effective":  85.0,
    "amp_pct_requested":   5.0,
    "sigma_tA":           -0.227,
    "R2":                  0.978,
    "tfold":              32.6,
    "label":              "Hartmann-resolved (warm-state)",
}



# Style

plt.rcParams.update({
    "axes.linewidth":      0.8,
    "axes.labelsize":     11,
    "xtick.labelsize":    10,
    "ytick.labelsize":    10,
    "xtick.direction":    "in",
    "ytick.direction":    "in",
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "legend.fontsize":    9,
    "legend.frameon":     False,
    "font.family":        "sans-serif",
    "mathtext.fontset":   "stix",
})


def main():
    outdir = Path("results")
    outdir.mkdir(exist_ok=True)


    fig, axs = plt.subplots(1, 2, figsize=(11, 4.5))

    ax = axs[0]

    # Debug mesh: connected scatter with R^2 colormap
    sc = ax.scatter(
        DEBUG["amp_pct"], DEBUG["sigma_tA"],
        c=DEBUG["R2"], s=130, cmap="viridis",
        vmin=0.7, vmax=1.0, edgecolors="k", linewidths=0.7,
        zorder=3, label="Debug mesh (2.30 M cells, same seed)",
    )
    ax.plot(DEBUG["amp_pct"], DEBUG["sigma_tA"],
            color="0.4", lw=0.9, zorder=2)

    # Hartmann-resolved single point (at effective amplitude)
    ax.scatter(
        [HRES["amp_pct_effective"]], [HRES["sigma_tA"]],
        marker="*", s=400, color="#cc3333",
        edgecolors="k", linewidths=0.9, zorder=4,
        label=f"{HRES['label']}\n(effective amplitude)",
    )

    # Annotate the Hartmann-resolved point
    ax.annotate(
        f"{HRES['sigma_tA']:+.3f}\n($R^2={HRES['R2']:.2f}$)",
        (HRES["amp_pct_effective"], HRES["sigma_tA"]),
        textcoords="offset points", xytext=(-95, 5),
        fontsize=9, color="#cc3333",
    )

    # Reference lines
    ax.axhline(0.0, color="k", lw=0.6)
    ax.axhline(-1.0, color="C3", lw=0.7, ls=":",
               label=r"$\sigma\tau_A = -1$ (Hartmann body force)")

    ax.set_xlabel(r"$|U'|_{\rm RMS}/|U|_{\rm lam}$ at probes (%)")
    ax.set_ylabel(r"$\sigma\,\tau_A$")
    ax.set_title(r"(a)  Mesh-sensitivity comparison: $\sigma\tau_A$ vs amplitude")
    ax.set_xlim(0, 100)
    ax.set_ylim(-1.05, 0.15)
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(loc="lower left", fontsize=9)

    cb = plt.colorbar(sc, ax=ax, fraction=0.045, pad=0.02)
    cb.set_label(r"$R^{\,2}$", fontsize=10)
    cb.ax.tick_params(labelsize=9)


    ax = axs[1]
    ax.plot(DEBUG["amp_pct"], DEBUG["R2"], "o-",
            color="#0d4f8b", lw=1.5, ms=8,
            label="Debug mesh")
    ax.scatter(
        [HRES["amp_pct_effective"]], [HRES["R2"]],
        marker="*", s=400, color="#cc3333",
        edgecolors="k", linewidths=0.9, zorder=4,
        label="Hartmann-resolved",
    )
    ax.set_xlabel(r"$|U'|_{\rm RMS}/|U|_{\rm lam}$ at probes (%)")
    ax.set_ylabel(r"$R^{\,2}$ of $\log\!{\rm RMS}\,U_z$ vs $t$")
    ax.set_title("(b)  Single-exponential fit quality vs amplitude")
    ax.set_xlim(0, 100)
    ax.set_ylim(0.75, 1.02)
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(loc="lower right", fontsize=9)

    # ----- Save -----
    fig.tight_layout()
    out_png = outdir / "mesh_comparison.png"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_png}")

  
    summary = outdir / "mesh_comparison_summary.md"
    with open(summary, "w") as f:
        f.write("# Mesh-sensitivity comparison: debug vs Hartmann-resolved\n\n")
        f.write("Caption"     )
        f.write("## Data table\n\n")
        f.write("| mesh | amplitude | $\\sigma\\tau_A$ | $R^{2}$ | "
                "e-fold (ms) |\n")
        f.write("|------|-----------|-----------------|---------|"
                "-------------|\n")
        for i in range(len(DEBUG["amp_pct"])):
            f.write(
                f"| debug | {DEBUG['amp_pct'][i]:.0f}\\,\\% (req.\\,=\\,probe) "
                f"| ${DEBUG['sigma_tA'][i]:+.3f}$ | "
                f"{DEBUG['R2'][i]:.2f} | {DEBUG['tfold'][i]:.1f} |\n"
            )
        f.write(
            f"| Hartmann-resolved | {HRES['amp_pct_effective']:.0f}\\,\\% "
            f"(effective; requested {HRES['amp_pct_requested']:.0f}\\,\\%) "
            f"| ${HRES['sigma_tA']:+.3f}$ | "
            f"{HRES['R2']:.2f} | {HRES['tfold']:.1f} |\n"
        )
    print(f"Wrote {summary}")


if __name__ == "__main__":
    main()
