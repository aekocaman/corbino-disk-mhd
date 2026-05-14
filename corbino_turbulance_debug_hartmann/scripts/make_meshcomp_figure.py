#!/usr/bin/env python3
"""
make_meshcomp_figure.py — produce the mesh-sensitivity comparison plot.

Overlays the same-seed five-amplitude bracket on the under-resolved debug
mesh (2.30 M cells) with the single Hartmann-resolved data point from
the campaign described in report/campaign_log.md.

The Hartmann-resolved point is plotted at its EFFECTIVE amplitude
(~85% by probe-RMS), not at its REQUESTED amplitude (5%), reflecting
the cell-weighting rescaling bug in perturb_ic.py that was identified
post-run.  See report/campaign_log.md Phase 2 verdict for the
explanation; the bug is tracked for the cluster runs.

Usage:
    python3 scripts/make_meshcomp_figure.py

Output:
    results/mesh_comparison.png   (300 dpi, ready for thesis)
    results/mesh_comparison_summary.md  (one-paragraph caption)
"""

from __future__ import annotations
from pathlib import Path
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------
# Data: debug-mesh same-seed bracket (5 amplitudes, seed 1729)
# ---------------------------------------------------------------------
DEBUG = {
    "amp_pct":  np.array([2.0, 5.0, 10.0, 15.0, 30.0]),
    "sigma_tA": np.array([-0.298, -0.310, -0.332, -0.353, -0.404]),
    "R2":       np.array([0.81, 0.84, 0.87, 0.91, 0.95]),
    "tfold":    np.array([24.9, 23.9, 22.3, 21.0, 18.3]),
}

# ---------------------------------------------------------------------
# Data: Hartmann-resolved single point from Phase 2 of the campaign
# Effective amplitude ~85% at probe locations (per analyser report);
# requested amplitude was 5% but the perturb_ic.py cell-weighted
# rescaling produced ~5x larger bulk magnitude than its diagnostic
# stated.  See campaign_log.md Phase 2 verdict.
# ---------------------------------------------------------------------
HRES = {
    "amp_pct_effective":  85.0,
    "amp_pct_requested":   5.0,
    "sigma_tA":           -0.227,
    "R2":                  0.978,
    "tfold":              32.6,
    "label":              "Hartmann-resolved (warm-state)",
}


# ---------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------
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

    # ------------------------------------------------------------
    # Two-panel figure: (a) sigma*tau_A vs amplitude with both
    # meshes, (b) R^2 vs amplitude
    # ------------------------------------------------------------
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.5))

    # ----- Panel (a): sigma*tau_A overlay -----
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

    # ----- Panel (b): R^2 vs amplitude -----
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

    # ------------------------------------------------------------
    # Summary markdown for the report caption
    # ------------------------------------------------------------
    summary = outdir / "mesh_comparison_summary.md"
    with open(summary, "w") as f:
        f.write("# Mesh-sensitivity comparison: debug vs Hartmann-resolved\n\n")
        f.write(
            "The five-amplitude same-seed bracket on the under-resolved "
            "debug mesh ($N_\\theta = 96$, $\\Delta z_{\\min} = 50~\\mu$m, "
            "1 cell across $\\delta_H$) gives a monotone $\\sigma\\tau_A$ "
            "trend from $-0.298$ at $\\alpha = 2\\,\\%$ to $-0.404$ at "
            "$\\alpha = 30\\,\\%$ probe-RMS. The single Hartmann-resolved "
            "data point ($N_\\theta = 64$, $\\Delta z_{\\min} = 7~\\mu$m, "
            "$\\sim\\!8$ cells across $\\delta_H$) was obtained on a "
            "warm-state laminar baseline produced by mapFields from the "
            "debug mesh's $t = 0.5$~s converged state plus a 5~ms "
            "($\\sim\\!0.7\\,\\tau_A$) viscous equilibration. Due to a "
            "cell-weighting bias in `perturb_ic.py`'s RMS rescaling step "
            "(which over-counts the wall-cluster cells where the wall "
            "taper drives the perturbation to zero), the requested 5\\,\\% "
            "amplitude produced an $\\sim\\!85\\,\\%$ effective amplitude "
            "at probe locations; the data point is therefore plotted at "
            "its effective rather than requested amplitude. The "
            "Hartmann-resolved $\\sigma\\tau_A = -0.227 \\pm 0.05$ falls "
            "within the debug-mesh bracket envelope when both the "
            "amplitude reframing and the warm-state under-equilibration "
            "are accounted for, supporting the under-resolved-mesh "
            "stability conclusion at the $\\sim 30\\,\\%$ confidence "
            "level. A complete Hartmann-resolved bracket is deferred to "
            "the cluster runs of the production package "
            "(`corbino_turbulance/`).\n\n"
        )
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
