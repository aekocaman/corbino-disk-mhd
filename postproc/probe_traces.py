#!/usr/bin/env python3
"""
Concatenates the original-run and the restart-run
probe histories at three radial probes (r = 4, 6, 8 cm, mid-plane,
φ = 0). Produces four figures:

  probe_traces.png            U_θ(t) + single-exp fit
                             

  probe_components.png        three stacked panels: U_θ(t), U_r(t),
                              U_z(t). U_r and U_z should be at the
                              numerical noise floor — sets the
                              axisymmetry threshold for Step-2
                              instability detection

  probe_two_exp_fit.png       same data with two-exponential fit
                                U_θ(t) = U_∞
                                       − A_A · exp(−t/τ_A)
                                       − A_H · exp(−t/τ_H)
                              τ_A locked at 7.5 ms (analytical bulk
                              Alfvén time); U_∞, A_A, A_H, τ_H are
                              the four free parameters. Reports
                              fitted τ_H per probe vs predicted
                              τ_H = 0.71 s.

  probe_rUtheta_collapse.png  r · U_θ(t) for all three probes
                              overlaid on one axis. Tests whether
                              the 1/r law is established
                              dynamically (curves overlie from
                              t ≳ τ_A) or only asymptotically
                              (curves merge near τ_H).

Usage:
    python3 probe_traces.py [/path/to/case]
"""

import sys
import re
import numpy as np
from pathlib import Path

PHYS_TAU_A = 7.5e-3   # bulk Alfvén time, s — predicted, locked in two-exp fit
PHYS_TAU_H = 0.71     # Hartmann time, s    — predicted, target of two-exp fit

# probe positions (must match system/probes_function or controlDict)
PROBES = {
    0: (0.04,  0.0,    0.005),   # r = 4 cm
    1: (0.06,  0.0,    0.005),   # r = 6 cm
    2: (0.08,  0.0,    0.005),   # r = 8 cm
}


def parse_probe_file(path: Path):
    """Parse OpenFOAM 'U' probe file and return (t, U_array[N,3]) for
    each probe column. The OpenFOAM file format is:
        # Probe N at (...)
        ...
        # Time   probe0  probe1  ...
        t   (Ux Uy Uz)   (Ux Uy Uz)   ...
    """
    text = path.read_text()
    # Strip header comment lines starting with '#'
    data_lines = [ln for ln in text.splitlines()
                  if ln.strip() and not ln.strip().startswith("#")]

    # Each line: t (Ux Uy Uz) (Ux Uy Uz) ...
    pat = re.compile(r"\(([-\d\.eE+]+)\s+([-\d\.eE+]+)\s+([-\d\.eE+]+)\)")
    times, U_per_probe = [], None
    for ln in data_lines:
        # split on whitespace, but tuples must stay together
        # easier: extract tuples then prepend the time
        toks = ln.split(None, 1)
        t = float(toks[0])
        rest = toks[1] if len(toks) > 1 else ""
        triples = pat.findall(rest)
        if U_per_probe is None:
            U_per_probe = [[] for _ in triples]
        for i, m in enumerate(triples):
            U_per_probe[i].append([float(m[0]), float(m[1]), float(m[2])])
        times.append(t)
    return np.array(times), np.array(U_per_probe)  # shape (N_probes, N_t, 3)


def az_unit_at(p):
    x, y, _ = p
    r = np.hypot(x, y)
    return np.array([-y / r, x / r, 0.0])


def rad_unit_at(p):
    x, y, _ = p
    r = np.hypot(x, y)
    return np.array([x / r, y / r, 0.0])


# ẑ unit vector is just (0, 0, 1) — projecting onto z gives U_z directly.


def main(case_dir: Path):
    pp = case_dir / "postProcessing" / "probes"
    if not pp.exists():
        print(f"[err] no postProcessing/probes/ under {case_dir}")
        sys.exit(1)

    # Find every numeric subdir (run start times) and sort
    run_dirs = sorted(
        [d for d in pp.iterdir()
         if d.is_dir() and re.fullmatch(r"[\d\.eE+\-]+", d.name)],
        key=lambda d: float(d.name),
    )
    print(f"[load] probe runs: {[d.name for d in run_dirs]}")

    # Stitch U files in chronological order
    all_t, all_U = [], None
    for d in run_dirs:
        u_file = d / "U"
        if not u_file.exists():
            continue
        t, U = parse_probe_file(u_file)
        if all_U is None:
            all_t, all_U = list(t), [list(seg) for seg in U]
        else:
            # avoid duplicate first sample of restart
            mask = t > all_t[-1]
            all_t.extend(t[mask])
            for i, seg in enumerate(U):
                all_U[i].extend(seg[mask])

    t = np.array(all_t)
    U = np.array(all_U)  # (N_probes, N_t, 3)
    print(f"[load] {len(t)} samples, {U.shape[0]} probes, "
          f"t ∈ [{t.min():.3f}, {t.max():.3f}] s")

    # Project U onto (θ̂, r̂, ẑ) at each probe location.
    # U_z is just the z-component; r̂ depends on probe position.
    U_theta = np.zeros((U.shape[0], U.shape[1]))
    U_r     = np.zeros((U.shape[0], U.shape[1]))
    U_z     = np.zeros((U.shape[0], U.shape[1]))
    for i, pid in enumerate(sorted(PROBES)):
        that = az_unit_at(PROBES[pid])
        rhat = rad_unit_at(PROBES[pid])
        U_theta[i] = U[i] @ that
        U_r[i]     = U[i] @ rhat
        U_z[i]     = U[i, :, 2]

    # Quick numerical sanity print: noise floors of U_r, U_z relative to U_θ
    print()
    print("Late-time (t > 1.0 s) component magnitudes — noise-floor diagnostic:")
    late = t > 1.0
    if late.any():
        print(f"{'probe':>5} {'r [cm]':>8} "
              f"{'|U_θ|':>11} {'|U_r|':>11} {'|U_z|':>11} "
              f"{'|U_r|/|U_θ|':>13} {'|U_z|/|U_θ|':>13}")
        print("-" * 80)
        for i, pid in enumerate(sorted(PROBES)):
            r_cm = np.hypot(*PROBES[pid][:2]) * 100
            mU_th = float(np.mean(np.abs(U_theta[i, late])))
            mU_r  = float(np.mean(np.abs(U_r[i, late])))
            mU_z  = float(np.mean(np.abs(U_z[i, late])))
            print(f"{pid:>5d} {r_cm:>8.1f} "
                  f"{mU_th:>11.3e} {mU_r:>11.3e} {mU_z:>11.3e} "
                  f"{mU_r/mU_th:>13.3e} {mU_z/mU_th:>13.3e}")

    # ------------------------------------------------------------------- #
    # Fit each probe's late-time tail to U_inf - A · exp(-t/τ)
    # ------------------------------------------------------------------- #
    from scipy.optimize import curve_fit
    def model(tt, U_inf, A, tau):
        return U_inf - A * np.exp(-tt / tau)

    print()
    print(f"{'probe':>5} {'r [cm]':>8} {'U_∞ fit [m/s]':>15} "
          f"{'τ fit [s]':>12} {'τ/τ_H_pred':>12}")
    print("-" * 60)

    fits = []
    fit_start = 0.10  # exclude very early transient (≈ 1 τ_A)
    for i, pid in enumerate(sorted(PROBES)):
        mask = t >= fit_start
        try:
            popt, _ = curve_fit(
                model, t[mask], U_theta[i, mask],
                p0=(U_theta[i, -1], 0.1, PHYS_TAU_H),
                maxfev=5000,
            )
            U_inf, A, tau = popt
            r_cm = np.hypot(*PROBES[pid][:2]) * 100
            print(f"{pid:>5d} {r_cm:>8.1f} {U_inf:+15.5f} "
                  f"{tau:>12.4f} {tau/PHYS_TAU_H:>12.3f}")
            fits.append((pid, U_inf, A, tau))
        except Exception as exc:
            print(f"[fit] probe {pid}: failed ({exc})")
            fits.append((pid, None, None, None))

    # ------------------------------------------------------------------- #
    # Two-exponential fit:  U_θ(t) = U_∞ − A_A·exp(−t/τ_A) − A_H·exp(−t/τ_H)
    # τ_A locked at PHYS_TAU_A = 7.5 ms (analytical bulk Alfvén time);
    # U_∞, A_A, A_H, τ_H are the four free parameters per probe.
    # ------------------------------------------------------------------- #
    def two_exp(tt, U_inf, A_A, A_H, tau_H):
        return (U_inf
                - A_A * np.exp(-tt / PHYS_TAU_A)
                - A_H * np.exp(-tt / tau_H))

    print()
    print("Two-exponential fit  (τ_A locked at 7.5 ms):")
    print(f"{'probe':>5} {'r [cm]':>8} {'U_∞':>11} "
          f"{'A_A':>11} {'A_H':>11} "
          f"{'τ_H fit [s]':>13} {'τ_H/τ_H_pred':>14}")
    print("-" * 80)

    fits2 = []
    fit_start2 = 0.005  # include the inertial phase so A_A is constrained
    for i, pid in enumerate(sorted(PROBES)):
        mask = t >= fit_start2
        try:
            popt, _ = curve_fit(
                two_exp, t[mask], U_theta[i, mask],
                p0=(U_theta[i, -1], 0.05, 0.10, PHYS_TAU_H),
                maxfev=10000,
            )
            U_inf, A_A, A_H, tau_H = popt
            r_cm = np.hypot(*PROBES[pid][:2]) * 100
            print(f"{pid:>5d} {r_cm:>8.1f} {U_inf:+11.5f} "
                  f"{A_A:+11.4e} {A_H:+11.4e} "
                  f"{tau_H:>13.4f} {tau_H/PHYS_TAU_H:>14.3f}")
            fits2.append((pid, U_inf, A_A, A_H, tau_H))
        except Exception as exc:
            print(f"[fit2] probe {pid}: failed ({exc})")
            fits2.append((pid, None, None, None, None))

    # ------------------------------------------------------------------- #
    # Plot
    # ------------------------------------------------------------------- #
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        COLORS = {0: "#d62728", 1: "#1f77b4", 2: "#2ca02c"}
        fig, ax = plt.subplots(figsize=(10, 5))

        # Build legend labels combining "r = X cm" with the per-probe
        # single-exp fit τ values — folds the upper-left annotation
        # box into the legend so we have a single info source on the
        # right.
        fit_by_pid = {p: tau for p, _, _, tau in fits}
        legend_labels = {}
        for i, pid in enumerate(sorted(PROBES)):
            r_cm = np.hypot(*PROBES[pid][:2]) * 100
            tau_i = fit_by_pid.get(pid)
            if tau_i is not None:
                legend_labels[pid] = (
                    f"r = {r_cm:.0f} cm    "
                    rf"$\tau$ = {tau_i:.3f} s    "
                    rf"($\tau/\tau_H$ = {tau_i/PHYS_TAU_H:.2f})"
                )
            else:
                legend_labels[pid] = f"r = {r_cm:.0f} cm"

        for i, pid in enumerate(sorted(PROBES)):
            ax.plot(t, U_theta[i], "-", color=COLORS[pid], lw=1.2,
                    label=legend_labels[pid])

            if fits[i][1] is not None:
                U_inf, A, tau = fits[i][1:]
                tt = np.linspace(fit_start, t.max(), 400)
                ax.plot(tt, model(tt, U_inf, A, tau),
                        "--", color=COLORS[pid], lw=1.0, alpha=0.8)
                ax.axhline(U_inf, color=COLORS[pid], ls=":", alpha=0.4)

        ax.set_xlabel("t  [s]")
        ax.set_ylabel(r"$U_\theta$  [m/s]")
        ax.set_title(rf"Probe Trace ($\tau_H$ = {PHYS_TAU_H:.2f} s)")
        ax.legend(loc="upper right", framealpha=0.9, ncol=1, fontsize=9)
        ax.grid(alpha=0.3)

        fig.tight_layout()
        out_png = case_dir / "probe_traces.png"
        fig.savefig(out_png, dpi=140)
        print(f"\n[save] {out_png}")
        plt.close(fig)

        # --------------------------------------------------------------- #
        # Companion 1 — three-component time series (axisymmetry diag.)
        # --------------------------------------------------------------- #
        fig2, axs = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
        labels_y = (r"$U_\theta$  [m/s]", r"$U_r$  [m/s]", r"$U_z$  [m/s]")
        data_per_panel = (U_theta, U_r, U_z)

        # Pre-compute late-time noise-floor ratios so we can fold them
        # into the U_θ-panel legend (one entry per probe carries r AND
        # the |U_r|/|U_θ|, |U_z|/|U_θ| ratios — no separate text box).
        late = t > 1.0
        legend_labels = {}
        for i, pid in enumerate(sorted(PROBES)):
            r_cm = np.hypot(*PROBES[pid][:2]) * 100
            if late.any():
                mU_th = float(np.mean(np.abs(U_theta[i, late])))
                mU_r  = float(np.mean(np.abs(U_r[i, late])))
                mU_z  = float(np.mean(np.abs(U_z[i, late])))
                legend_labels[pid] = (
                    f"r = {r_cm:.0f} cm    "
                    rf"$|U_r|/|U_\theta|$ = {mU_r/mU_th:.1e}    "
                    rf"$|U_z|/|U_\theta|$ = {mU_z/mU_th:.1e}"
                )
            else:
                legend_labels[pid] = f"r = {r_cm:.0f} cm"

        for j, (axp, ylabel, data) in enumerate(
                zip(axs, labels_y, data_per_panel)):
            for i, pid in enumerate(sorted(PROBES)):
                # Only attach legend labels in the first (U_θ) panel.
                lbl = legend_labels[pid] if j == 0 else None
                axp.plot(t, data[i], "-", color=COLORS[pid], lw=1.0,
                         label=lbl)
            axp.set_ylabel(ylabel)
            axp.grid(alpha=0.3)
            axp.axhline(0.0, color="k", lw=0.5, alpha=0.4)

        # Single legend in the U_θ panel — colour ↔ r mapping plus the
        # noise-floor ratios that used to live in a separate text box.
        axs[0].legend(
            loc="upper right", framealpha=0.9, ncol=1, fontsize=9,
        )
        axs[0].set_title("Velocity Components at Probes")
        axs[2].set_xlabel("t  [s]")
        fig2.tight_layout()
        out_png2 = case_dir / "probe_components.png"
        fig2.savefig(out_png2, dpi=140)
        print(f"[save] {out_png2}")
        plt.close(fig2)

        # --------------------------------------------------------------- #
        # Companion 2 — two-exponential fit overlaid on U_θ(t)
        # --------------------------------------------------------------- #
        fig3, ax3 = plt.subplots(figsize=(10, 5))
        for i, pid in enumerate(sorted(PROBES)):
            r_cm = np.hypot(*PROBES[pid][:2]) * 100
            ax3.plot(t, U_theta[i], "-", color=COLORS[pid], lw=1.2,
                     label=f"r = {r_cm:.0f} cm")
            if fits2[i][1] is not None:
                _, U_inf, A_A, A_H, tau_H = fits2[i]
                tt = np.linspace(fit_start2, t.max(), 800)
                ax3.plot(tt, two_exp(tt, U_inf, A_A, A_H, tau_H),
                         "--", color=COLORS[pid], lw=1.0, alpha=0.85)
                ax3.axhline(U_inf, color=COLORS[pid], ls=":", alpha=0.4)

        ax3.set_xlabel("t  [s]")
        ax3.set_ylabel(r"$U_\theta$  [m/s]")
        ax3.set_title(
            r"Two-exponential fit: $U_\infty - A_A\,e^{-t/\tau_A} "
            r"- A_H\,e^{-t/\tau_H}$,  "
            rf"$\tau_A$ = {PHYS_TAU_A*1e3:.1f} ms locked",
            pad=28,
        )
        ax3.legend(
            loc="lower center", bbox_to_anchor=(0.5, 1.02),
            ncol=3, frameon=False,
        )
        ax3.grid(alpha=0.3)
        txt2 = "\n".join(
            f"r = {np.hypot(*PROBES[p][:2])*100:.0f} cm:  "
            rf"$\tau_H$ = {tau_H:.3f} s   "
            rf"($\tau_H/\tau_H^{{\,\mathrm{{pred}}}}$ = "
            f"{tau_H/PHYS_TAU_H:.2f})"
            for p, _, _, _, tau_H in fits2 if tau_H is not None
        )
        ax3.text(0.02, 0.05, txt2, transform=ax3.transAxes,
                 fontsize=9, va="bottom", family="monospace",
                 bbox=dict(facecolor="white", alpha=0.85, edgecolor="0.7"))
        fig3.tight_layout()
        out_png3 = case_dir / "probe_two_exp_fit.png"
        fig3.savefig(out_png3, dpi=140)
        print(f"[save] {out_png3}")
        plt.close(fig3)

        # --------------------------------------------------------------- #
        # Companion 3 — r·U_θ(t) collapse (dynamic 1/r law)
        # --------------------------------------------------------------- #
        fig4, ax4 = plt.subplots(figsize=(10, 5))
        for i, pid in enumerate(sorted(PROBES)):
            r_cm = np.hypot(*PROBES[pid][:2]) * 100
            r_m  = np.hypot(*PROBES[pid][:2])
            ax4.plot(t, r_m * U_theta[i], "-", color=COLORS[pid], lw=1.2,
                     label=f"r = {r_cm:.0f} cm")
        # Late-time mean as reference
        late = t > 1.0
        if late.any():
            rUtheta_means = []
            for i, pid in enumerate(sorted(PROBES)):
                r_m = np.hypot(*PROBES[pid][:2])
                rUtheta_means.append(float(np.mean(r_m * U_theta[i, late])))
            mean_ref = float(np.mean(rUtheta_means))
            spread = (max(rUtheta_means) - min(rUtheta_means)) / abs(mean_ref)
            ax4.axhline(mean_ref, color="k", ls=":", alpha=0.5,
                        label=f"⟨r·U_θ⟩ at t > 1 s = {mean_ref:+.3e}")
            ax4.text(0.02, 0.05,
                     f"3-probe spread of r·U_θ at t > 1 s: "
                     f"{100*spread:.2f} %",
                     transform=ax4.transAxes, fontsize=9,
                     va="bottom", family="monospace",
                     bbox=dict(facecolor="white", alpha=0.85,
                               edgecolor="0.7"))
        ax4.set_xlabel("t  [s]")
        ax4.set_ylabel(r"$r \cdot U_\theta$  [m$^2$/s]")
        ax4.set_title(
            r"Dynamic 1/$r$ law — $r \cdot U_\theta(t)$ at three radii",
            pad=28,
        )
        ax4.legend(
            loc="lower center", bbox_to_anchor=(0.5, 1.02),
            ncol=4, frameon=False, fontsize=9,
        )
        ax4.grid(alpha=0.3)
        fig4.tight_layout()
        out_png4 = case_dir / "probe_rUtheta_collapse.png"
        fig4.savefig(out_png4, dpi=140)
        print(f"[save] {out_png4}")
        plt.close(fig4)

    except Exception as exc:
        print(f"[plot] skipped ({exc})")


if __name__ == "__main__":
    case = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd()
    main(case.resolve())
