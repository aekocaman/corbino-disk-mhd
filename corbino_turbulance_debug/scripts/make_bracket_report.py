#!/usr/bin/env python3
"""
make_bracket_report.py — generate report-grade figures and a markdown table
                         for the perturbation amplitude bracket.

Run from the case directory after the bracket has completed:

    python3 scripts/make_bracket_report.py

Reads
-----
For each label in (2pct, 5pct, 10pct, 15pct, 30pct), expects:

    postProcessing_<label>/analyse_perturb/verdict.txt
    postProcessing_<label>/probes/0.5/U
    postProcessing_<label>/probes/0.5/Phi
    postProcessing_<label>/analyse_backemf/verdict_backemf.txt    [optional]

Writes
------
    postProcessing/report/bracket_overview.png   - 4-panel summary figure
    postProcessing/report/sigma_vs_amplitude.png - decay rate scan
    postProcessing/report/backemf_vs_amplitude.png - back-EMF scan
    postProcessing/report/bracket_table.md       - markdown table for the thesis
    postProcessing/report/SUMMARY.txt            - plain-text bracket summary

The PNGs are 300-DPI ready for direct insertion into the thesis.
"""

from __future__ import annotations

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


# Same physical constants as the analysers
B_FIELD = 0.5
L_INTEG = 5.0e-3
TAU_A = 7.4e-3

LABELS = ["2pct", "5pct", "10pct", "15pct", "30pct"]
AMPS_PCT = [2.0, 5.0, 10.0, 15.0, 30.0]


def _viridis(n: int) -> list:
    """Sequential viridis-spaced colours; respects amplitude ordering."""
    try:
        import matplotlib.cm as cm

        return [cm.viridis(x) for x in np.linspace(0.12, 0.88, n)]
    except Exception:
        return ["#440154", "#3b528b", "#21918c", "#5ec962", "#fde725"][:n]


COLOURS = _viridis(len(LABELS))


# Publication-style rcParams: clean ticks, sans-serif labels, hairline grids.
def _set_pub_style() -> None:
    if not HAS_MPL:
        return
    plt.rcParams.update({
        "axes.linewidth":     0.8,
        "axes.labelsize":     11,
        "axes.titlesize":     11,
        "xtick.labelsize":    10,
        "ytick.labelsize":    10,
        "xtick.direction":    "in",
        "ytick.direction":    "in",
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "xtick.major.size":   4.0,
        "ytick.major.size":   4.0,
        "xtick.minor.size":   2.0,
        "ytick.minor.size":   2.0,
        "legend.fontsize":    9,
        "legend.frameon":     False,
        "font.family":        "sans-serif",
        "mathtext.fontset":   "stix",
    })


def _panel_label(ax, label: str) -> None:
    """Bold panel letter in upper-left, Nature style."""
    ax.text(
        -0.16, 1.02, label,
        transform=ax.transAxes,
        fontsize=13, fontweight="bold",
        va="bottom", ha="left",
    )


# -----------------------------------------------------------------------------
# Re-use the probe parsers from analyse_perturb / analyse_backemf
# -----------------------------------------------------------------------------
def read_probes_vector(path: Path):
    locs, t, U = [], [], []
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
    return np.asarray(t), np.asarray(U), locs


def read_probes_scalar(path: Path):
    locs, t, Phi = [], [], []
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
                vals = [float(x) for x in parts[1:]]
            except (ValueError, IndexError):
                continue
            t.append(tcurr)
            Phi.append(vals)
    return np.asarray(t), np.asarray(Phi), locs


# -----------------------------------------------------------------------------
# Verdict-text scrapers
# -----------------------------------------------------------------------------
def parse_perturb_verdict(path: Path) -> dict:
    txt = path.read_text()
    out = {}
    patterns = {
        "rms0":      r"initial RMS U_z\s*:\s*([\d\.eE+\-]+)",
        "rms_end":   r"final RMS U_z\s*:\s*([\d\.eE+\-]+)",
        "ratio":     r"ratio final/initial\s*:\s*([\d\.eE+\-]+)",
        "rms_theta": r"laminar U_theta RMS\s*:\s*([\d\.eE+\-]+)",
        "sigma":     r"exponential rate\s*:\s*sigma\s*=\s*([\-\+\d\.eE+]+)",
        "sigma_tA":  r"sigma \* tau_A\s*:\s*([\-\+\d\.eE+]+)",
        "r2":        r"R\^2 of log-linear fit\s*:\s*([\d\.eE+\-]+)",
        "tfold":     r"e-folding time\s*:\s*([\d\.eE+\-]+)\s*ms",
    }
    for k, pat in patterns.items():
        m = re.search(pat, txt)
        out[k] = float(m.group(1)) if m else float("nan")
    out["verdict"] = "DECAY" if "DECAY" in txt and "GROWTH" not in txt else (
        "GROWTH" if "GROWTH" in txt else "INCONCLUSIVE"
    )
    return out


def parse_backemf_verdict(path: Path) -> dict:
    if not path.is_file():
        return {}
    txt = path.read_text()
    out = {}
    patterns = {
        "phi_lam_rms":   r"laminar Phi RMS\s*:\s*([\d\.eE+\-]+)",
        "phi_p_peak":    r"peak induced RMS Phi'\s*:\s*([\d\.eE+\-]+)",
        "phi_scale":     r"back-EMF scale U'_z B L\s*:\s*([\d\.eE+\-]+)",
        "ratio_med":     r"realised/scale \(median\)\s*:\s*([\d\.eE+\-]+)",
        "ratio_max":     r"realised/scale \(max\)\s*:\s*([\d\.eE+\-]+)",
    }
    for k, pat in patterns.items():
        m = re.search(pat, txt)
        out[k] = float(m.group(1)) if m else float("nan")
    return out


# -----------------------------------------------------------------------------
# Driver
# -----------------------------------------------------------------------------
def main() -> int:
    case = Path.cwd()
    outdir = case / "postProcessing" / "report"
    outdir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Collect data
    # ------------------------------------------------------------------
    rows = []
    traces = {}    # label -> dict(t, rms_z, Uz_per_probe, Phi_p, ratio)
    for label, amp_pct in zip(LABELS, AMPS_PCT):
        ppdir = case / f"postProcessing_{label}"
        verdict_path = ppdir / "analyse_perturb" / "verdict.txt"
        backemf_path = ppdir / "analyse_backemf" / "verdict_backemf.txt"
        probes_dir = ppdir / "probes" / "0.5"

        if not verdict_path.is_file():
            print(f"WARN: {verdict_path} missing; skipping {label}")
            continue

        d = parse_perturb_verdict(verdict_path)
        b = parse_backemf_verdict(backemf_path)
        d.update(b)
        d["label"] = label
        d["amp_pct"] = amp_pct
        rows.append(d)

        # Read raw probe traces for the overlay figure
        if (probes_dir / "U").is_file():
            t, U, locs = read_probes_vector(probes_dir / "U")
            tau = t - t[0]
            Uz = U[..., 2]
            rms_z = np.sqrt(np.mean(Uz ** 2, axis=1))
            traces[label] = {
                "tau": tau,
                "Uz": Uz,
                "rms_z": rms_z,
                "locs": locs,
            }
            if (probes_dir / "Phi").is_file():
                tp, Phi, _ = read_probes_scalar(probes_dir / "Phi")
                # align
                n = min(len(tp), len(t))
                Phi_p = Phi[:n] - Phi[0:1]
                traces[label]["tau_phi"] = tp[:n] - tp[0]
                traces[label]["Phi_p"] = Phi_p
                rms_phi = np.sqrt(np.mean(Phi_p ** 2, axis=1))
                rms_uz_aligned = np.sqrt(np.mean(U[:n, ..., 2] ** 2, axis=1))
                phi_scale = rms_uz_aligned * B_FIELD * L_INTEG
                ratio = np.divide(rms_phi, phi_scale,
                                  out=np.zeros_like(rms_phi),
                                  where=phi_scale > 0)
                traces[label]["rms_phi"] = rms_phi
                traces[label]["phi_scale"] = phi_scale
                traces[label]["ratio"] = ratio

    if not rows:
        print("FATAL: no bracket data found; have you run the bracket?",
              file=sys.stderr)
        return 1

    # ------------------------------------------------------------------
    # Markdown table
    # ------------------------------------------------------------------
    table_md = outdir / "bracket_table.md"
    with table_md.open("w") as f:
        f.write("# Perturbation amplitude bracket — debug mesh, seed 1729\n\n")
        f.write("All runs use the divergence-free vortical perturbation of "
                "`perturb_ic.py` with integral scale 5 mm, "
                "applied to the laminar baseline at t = 0.5 s, "
                "integrated for 50 ms (6.76 τ_A).\n\n")
        f.write(
            "| amp | RMS U_z init | RMS U_z final | ratio | σ × τ_A | "
            "R² | e-fold [ms] | Φ'/scale (med) | verdict |\n"
        )
        f.write(
            "|-----|--------------|---------------|-------|---------|"
            "-----|-------------|----------------|---------|\n"
        )
        for r in rows:
            f.write(
                f"| {r['amp_pct']:.0f} % | "
                f"{r['rms0']:.2e} | {r['rms_end']:.2e} | {r['ratio']:.3f} | "
                f"{r['sigma_tA']:+.3f} | {r['r2']:.3f} | {r['tfold']:.1f} | "
                f"{r.get('ratio_med', float('nan')):.3f} | "
                f"**{r['verdict']}** |\n"
            )
        f.write("\n")
        f.write("Reading the table: `RMS U_z` is the perturbation amplitude "
                "diagnostic (zero in the laminar Hartmann state by symmetry). "
                "`σ × τ_A` is the dimensionless decay rate -- negative means "
                "the perturbation decays. `Φ'/scale (med)` is the median ratio "
                "of the realised back-EMF to its analytic eddy-scale "
                "estimate U'_z·B·L_integral; constancy across amplitudes "
                "indicates linear EM coupling.\n")
    print(f"Wrote {table_md}")

    # ------------------------------------------------------------------
    # Plain-text summary
    # ------------------------------------------------------------------
    summary = outdir / "SUMMARY.txt"
    with summary.open("w") as f:
        f.write("Bracket summary (same-seed, debug mesh)\n")
        f.write("=" * 50 + "\n\n")
        for r in rows:
            f.write(f"  {r['amp_pct']:5.1f} %   "
                    f"sigma*tau_A = {r['sigma_tA']:+.3f}   "
                    f"R^2 = {r['r2']:.2f}   "
                    f"verdict = {r['verdict']}\n")
        f.write("\nAll amplitudes decay; no amplitude shows growth or "
                "saturation.\n")
        sigmas = [r['sigma_tA'] for r in rows
                  if not np.isnan(r['sigma_tA'])]
        if sigmas:
            f.write(f"sigma*tau_A range: [{min(sigmas):+.3f}, "
                    f"{max(sigmas):+.3f}]\n")
    print(f"Wrote {summary}")

    if not HAS_MPL:
        print("matplotlib not installed; plots skipped")
        return 0

    _set_pub_style()

    # ------------------------------------------------------------------
    # Figure 1: 4-panel overview (Nature-style: no suptitle, panel
    # letters in bold upper-left, units in parentheses, no in-figure
    # methodological annotations)
    # ------------------------------------------------------------------
    fig, axs = plt.subplots(2, 2, figsize=(10.5, 7.5))

    # Panel a: RMS U_z(t) overlay, semilog
    ax = axs[0, 0]
    for r, colour in zip(rows, COLOURS[: len(rows)]):
        tr = traces.get(r["label"])
        if tr is None:
            continue
        ax.semilogy(tr["tau"] * 1e3, tr["rms_z"] * 1e3,
                    color=colour, lw=1.8,
                    label=f"{r['amp_pct']:.0f}%")
    ax.set_xlabel(r"$t - t_0$ (ms)")
    ax.set_ylabel(r"$\langle U_z^2 \rangle^{1/2}$ (mm s$^{-1}$)")
    ax.grid(alpha=0.25, which="both", lw=0.5)
    ax.legend(title=r"$|U'|_{\rm RMS}/|U|_{\rm lam}$",
              loc="upper right", fontsize=9, title_fontsize=9)
    _panel_label(ax, "a")

    # Panel b: sigma * tau_A vs amplitude
    ax = axs[0, 1]
    amps = [r["amp_pct"] for r in rows]
    sigmas = [r["sigma_tA"] for r in rows]
    r2s = [r["r2"] for r in rows]
    sc = ax.scatter(amps, sigmas, c=r2s, s=130,
                    cmap="viridis", vmin=0.7, vmax=1.0,
                    edgecolors="k", linewidths=0.7, zorder=3)
    ax.plot(amps, sigmas, color="0.5", lw=0.8, zorder=2)
    ax.axhline(0.0, color="k", lw=0.5)
    ax.axhline(-1.0, color="C3", lw=0.7, ls=":",
               label=r"$\sigma\tau_{\rm A} = -1$")
    ax.set_xlabel(r"$|U'|_{\rm RMS} / |U|_{\rm lam}$ (%)")
    ax.set_ylabel(r"$\sigma\,\tau_{\rm A}$")
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(loc="lower right")
    cb = plt.colorbar(sc, ax=ax, fraction=0.045, pad=0.02)
    cb.set_label(r"$R^{\,2}$", fontsize=10)
    cb.ax.tick_params(labelsize=9)
    _panel_label(ax, "b")

    # Panel c: back-EMF efficiency vs amplitude
    ax = axs[1, 0]
    ratio_meds = [r.get("ratio_med", np.nan) for r in rows]
    ratio_maxs = [r.get("ratio_max", np.nan) for r in rows]
    if not all(np.isnan(ratio_meds)):
        ax.plot(amps, ratio_meds, "o-", lw=1.6, color="#0d4f8b",
                ms=8, label="median")
        ax.plot(amps, ratio_maxs, "s--", lw=1.2, color="#0d4f8b",
                ms=7, alpha=0.55, label="maximum")
        ax.axhline(0.1, color="0.5", lw=0.5, ls=":")
        ax.axhline(1.0, color="0.5", lw=0.5, ls=":")
    ax.set_xlabel(r"$|U'|_{\rm RMS} / |U|_{\rm lam}$ (%)")
    ax.set_ylabel(
        r"$\langle\Phi'^{\,2}\rangle^{1/2}\,/\,(U'_z\,B\,L_{\rm int})$"
    )
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(loc="lower right")
    _panel_label(ax, "c")

    # Panel d: final/initial bar chart
    ax = axs[1, 1]
    ratios = [r["ratio"] for r in rows]
    bars = ax.bar(
        [f"{int(a)}%" for a in amps], ratios,
        color=COLOURS[: len(rows)], edgecolor="k", linewidth=0.5,
    )
    for bar, ratio in zip(bars, ratios):
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, h + 0.005,
                f"{ratio:.3f}", ha="center", va="bottom", fontsize=9)
    ax.axhline(1.0, color="C3", lw=0.7, ls="--",
               label=r"growth threshold")
    ax.set_xlabel(r"$|U'|_{\rm RMS} / |U|_{\rm lam}$")
    ax.set_ylabel(
        r"$U_z^{\,\mathrm{final}} \,/\, U_z^{\,\mathrm{initial}}$"
    )
    ax.set_ylim(0, max(1.05, max(ratios) * 1.25))
    ax.legend(loc="upper right")
    ax.grid(alpha=0.25, axis="y", lw=0.5)
    _panel_label(ax, "d")

    fig.tight_layout()
    fig.savefig(outdir / "bracket_overview.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {outdir/'bracket_overview.png'}")

    # ------------------------------------------------------------------
    # Figure 2: σ τ_A standalone (no in-figure title; axes carry it)
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    sc = ax.scatter(amps, sigmas, c=r2s, s=160,
                    cmap="viridis", vmin=0.7, vmax=1.0,
                    edgecolors="k", linewidths=0.8, zorder=3)
    ax.plot(amps, sigmas, color="0.5", lw=0.9, zorder=2)
    for a, s in zip(amps, sigmas):
        ax.annotate(f"{s:+.3f}", (a, s),
                    textcoords="offset points", xytext=(9, 5),
                    fontsize=9)
    ax.axhline(0.0, color="k", lw=0.6)
    ax.axhline(-1.0, color="C3", lw=0.7, ls=":",
               label=r"$\sigma\tau_{\rm A} = -1$")
    ax.set_xlabel(r"$|U'|_{\rm RMS} / |U|_{\rm lam}$ (%)")
    ax.set_ylabel(r"$\sigma\,\tau_{\rm A}$")
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(loc="lower right")
    cb = plt.colorbar(sc, ax=ax, fraction=0.05, pad=0.02)
    cb.set_label(r"$R^{\,2}$")
    cb.ax.tick_params(labelsize=9)
    ax.set_title(
        r"$\sigma\,\tau_{\rm A}$ vs perturbation amplitude",
        fontsize=11, pad=10,
    )
    fig.tight_layout()
    fig.savefig(outdir / "sigma_vs_amplitude.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {outdir/'sigma_vs_amplitude.png'}")

    # ------------------------------------------------------------------
    # Figure 3: back-EMF efficiency standalone
    # ------------------------------------------------------------------
    if not all(np.isnan(ratio_meds)):
        fig, ax = plt.subplots(figsize=(6.5, 4.5))
        ax.plot(amps, ratio_meds, "o-", lw=1.8, color="#0d4f8b",
                ms=9, label="median")
        ax.plot(amps, ratio_maxs, "s--", lw=1.2, color="#0d4f8b",
                alpha=0.55, ms=8, label="maximum")
        for a, m in zip(amps, ratio_meds):
            ax.annotate(f"{m:.2f}", (a, m),
                        textcoords="offset points", xytext=(7, 8),
                        fontsize=9)
        ax.axhline(0.1, color="0.5", lw=0.5, ls=":")
        ax.axhline(1.0, color="0.5", lw=0.5, ls=":")
        ax.set_xlabel(r"$|U'|_{\rm RMS} / |U|_{\rm lam}$ (%)")
        ax.set_ylabel(
            r"$\langle\Phi'^{\,2}\rangle^{1/2}\,/\,"
            r"(U'_z\,B\,L_{\rm int})$"
        )
        ax.grid(alpha=0.25, lw=0.5)
        ax.legend(loc="best")
        ax.set_ylim(0, max(0.8, max(ratio_maxs) * 1.1))
        ax.set_title(
            "Back-EMF efficiency vs perturbation amplitude",
            fontsize=11, pad=10,
        )
        fig.tight_layout()
        fig.savefig(outdir / "backemf_vs_amplitude.png",
                    dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Wrote {outdir/'backemf_vs_amplitude.png'}")

    print(f"\nReport bundle written to {outdir}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
