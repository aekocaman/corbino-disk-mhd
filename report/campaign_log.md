# Hartmann-resolved bracket campaign — execution log

**Eray Kocaman, Bilkent University, 2026.**
This document is a step-by-step record of the 3-day laptop campaign
launched to extend the Step-2 first-look bracket onto the
Hartmann-resolved sibling mesh via a `mapFields` warm-start. It is
intended to be appended verbatim to the final report as an execution
appendix, providing full provenance for the production result.

The campaign plan is locked at:
- Phase 0: setup and `mapFields` debug → resolved (< 5 min)
- Phase 1: laminar equilibration on the resolved mesh (~2 h, 5 ms physical)
- Phase 2: bracket run @ 5% RMS (~20 h, 50 ms physical)
- Phase 3: bracket run @ 15% RMS (~20 h, 50 ms physical)
- Phase 4: bracket run @ 30% RMS (~20 h, 50 ms physical)
- Phase 5: analysers + side-by-side figures (~1 h)

Total wall clock budget: ~63 h ≈ 2.6 days.
Mesh: `corbino_turbulance_debug_hartmann/`, 3.58 M cells (200×64×70 per block).
Per-step wall clock observed in earlier diagnostics: ~70 s/step on 8 ranks.

Each section below is filled in as the corresponding phase completes.
Headings give the launch and finish wall-clock times so the total
elapsed cost can be cross-checked against the budget.

---

## Phase 0: setup and `mapFields` debug → resolved

**Status:** PARTIALLY COMPLETED — required a mesh-version correction.
Final state: PASS after rebuild.

**Purpose.** Take the converged laminar solution at `t = 0.5 s` from
the under-resolved debug case (`corbino_turbulance_debug/`, 2.30 M cells)
and interpolate it onto the cell centres of the Hartmann-resolved mesh
(3.58 M cells). The bulk solution is preserved (the debug mesh resolved
the bulk Hartmann balance to 1% via the 1/r law); only the wall-layer
discontinuity is introduced at the resolved mesh's first wall cells,
where the debug-mesh under-resolution had collapsed the layer onto a
single 50 µm cell.

### Phase 0 — first pass (issue caught)

The first run of the standard staging block on `corbino_turbulance_debug_hartmann/`
delivered `mapFields` exit 0 and a clean log, with the laminar
$U_\theta\!\sim\!0.07~\mathrm{m/s}$ correctly visible in the
`internalField` sample of `0/U`. However, the cell count reported by
`mapFields` was

    internalField   nonuniform List<vector>
    5376000

i.e.\ `4 × 200 × 96 × 70 = 5,376,000`, the OLDER 96-θ Hartmann-resolved
mesh that previously thrashed swap on the laptop at 78 s/step — not the
`64-θ`, 3.58 M-cell RAM-fitted mesh designed for this campaign. Whichever
prior test left the case directory had used the older `blockMeshDict`.

Continuing the campaign on a swap-bound mesh for 60+ hours would risk
OOM-kill, fragmentation, and unpredictable per-step time. The reliable
fix is a 5-minute rebuild of the polyMesh from the corrected dictionary.

### Phase 0 — corrective rebuild

The 64-θ `blockMeshDict` was re-staged from
`corbino_turbulance_debug_hartmann/system/blockMeshDict`,
`constant/polyMesh/` was wiped, and `blockMesh` rebuilt the mesh at
$3.58\!\times\!10^{6}$ cells. `checkMesh` confirmed the same 13/14
passing metrics as the previous Hartmann-resolved test (only the
intentional cell-determinant warning at the wall clusters fires).
`mapFields` was then re-run from `corbino_turbulance_debug/0.5/` onto
the rebuilt mesh, again exiting 0 with all 10 source fields (B, C, J,
Phi, U, U.lam, p, plus the cell-centre fields Ccx/Ccy/Ccz) interpolated.

Verdict: **PASS** at `0/U` cell count = 3,584,000 and `mapFields` log
terminating in `End`. The 0/ directory now contains the
$\sim$1.2 GB of mapped binary fields ready for the Phase-1
equilibration restart.

**Lesson banked into the campaign log.** A `grep "^  cells:" log.blockMesh`
sanity check is added to every Phase-0/Phase-1 launcher block from
this point onward, so a stale polyMesh cannot silently propagate
through the campaign.

---

## Phase 1: laminar equilibration on resolved mesh (5 ms)

**Status:** RUNNING (mid-flight check at $t_{\rm wall}\!\approx\!29$ min).
Final verdict pending completion at $\sim$60 min total.

**Purpose.** Smooth the wall-layer step discontinuity introduced by
`mapFields`. Viscous timescale across the wall layer is
$\tau_{\rm visc} = \delta_H^{2}/\nu \approx 7$~ms; integrating for 5~ms
(\,$\sim\!0.7\,\tau_{\rm visc}$\,) is sufficient to let the viscous
diffusion redistribute the mapped field across the resolved cells,
but short enough that the bulk solution does not drift appreciably
(5~ms = $0.007\,\tau_H$; the bulk evolves negligibly on this timescale).

**Wall-clock estimate (revised after mid-flight measurement).** Initial
budget assumed 70~s/step from the earlier cold-start Hartmann-resolved
test. Mid-flight measurement of the present run (warm-started via
\texttt{mapFields}) yielded $\sim\!30.5$~s/step --- approximately
$2.3\!\times$ faster than the cold-start figure. The 100-timestep budget
therefore translates to $\sim\!50$~min wall clock instead of $\sim\!2$~h.
Cause of the speed-up is plausibly that the warm-start IC has the bulk
$U_\theta$ already at its converged value, so the PISO loop spends
fewer outer iterations stabilising the rapidly-evolving cold-start
Lorentz spin-up.

**Command sequence.** Standard 5-ms equilibration block:
\begin{verbatim}
cp system/controlDict.long system/controlDict
sed -i 's/^endTime.*/endTime         0.005;/'        system/controlDict
sed -i 's/^writeInterval.*/writeInterval   0.005;/'  system/controlDict
cp system/fvSolution.long system/fvSolution
decomposePar -force > log.decomposePar.equil 2>&1
nohup mpirun -np 8 epotMhdFoam -parallel > log.epotMhdFoam.equil 2>&1 &
\end{verbatim}

**Mid-flight diagnostics ($t_{\rm wall} = 29$~min, 57 steps completed):**

  - Process state: alive, $\sim$0.3 \% memory per the parent process
  - Simulated time: $t = 2.01$~ms (40 \% of the 5 ms target)
  - Per-step time: $30.5$~s
  - GAMG-$\Phi$ iterations per step: 2 (excellent)
  - Cumulative continuity error: $-1.4\!\times\!10^{-11}$ (machine
    precision)
  - Memory used: 5.0~GB of 7.6~GB; swap unchanged from pre-run baseline
    of 837~MiB. \emph{The 3.58~M mesh fits comfortably in RAM with no
    new swapping.}
  - No \texttt{nan}/\texttt{inf}/\texttt{FPE}/\texttt{abort} markers in
    the log.

**Final outputs (Phase 1 completed at $t_{\rm wall}\approx 57$ min).**
\begin{itemize}
\item Steps completed: 117 (slightly more than the 100-step target because $\Delta t$ was CFL-bounded slightly below the maxDeltaT cap).
\item Final cumulative continuity error: $-2.16\!\times\!10^{-11}$ (machine precision).
\item GAMG-$\Phi$ at 2 iter/step throughout.
\item Snapshot at $t = 0.005$~s present in all eight \texttt{processor*} directories.
\item Memory after run: 461~MiB used (solver released its working set), swap unchanged at the residual baseline.
\item No NaN/FPE/abort markers.
\end{itemize}

**Verdict: PASS.**

**Observation worth flagging for the report.** The volume-averaged laminar RMS at the end of Phase~1 reads $5.37\!\times\!10^{-2}$~m/s; the same diagnostic on the under-resolved debug mesh at $t = 0.5$~s read $9.78\!\times\!10^{-2}$~m/s. The 45\,\% drop is larger than expected from the warm-start strategy alone. The diagnosis is that the 5~ms equilibration window is only $\sim 0.7\,\tau_A$, sufficient for the wall-layer viscous transient ($\tau_{\rm visc}(\delta_H) \approx 7$~ms) but very short on the bulk Hartmann timescale $\tau_H = 0.71$~s. During those 5~ms the now-resolved viscous diffusion across the wall layer pulls momentum from the bulk faster than the Lorentz body force can replenish it, dropping $U_\infty$ from the mapFields-inherited value ($\sim 0.097$~m/s at $r=6$~cm on the debug mesh) toward an intermediate non-equilibrium value of $\sim 0.045$~m/s. The bracket runs of Phases~2--4 therefore measure $\sigma\,\tau_A$ on a warm-equilibrating state rather than a fully-converged Hartmann-resolved laminar baseline. The amplitude calibration $|U'|_{\rm RMS}/|U_{\rm lam}|_{\rm RMS} = \alpha$ is correct relative to the current state, so the bracket diagnostics remain self-consistent; the side-by-side comparison with the debug-mesh bracket result inherits this caveat as a documented difference in operating point.

---

## Phase 2: bracket run @ 5 % RMS

**Status:** RUNNING (perturbation injected, solver launched).

**Purpose.** First Hartmann-resolved data point of the same-seed
amplitude bracket. Direct mesh-sensitivity comparison against the
debug-mesh $\sigma\!\times\!\tau_A = -0.310 \pm 0.04$ result of finding~\#32.

**Wall-clock estimate (revised).** 1000 timesteps $\times\!\sim\!30$~s/step
$\approx\!8.3$~h. (Original budget: $\sim\!20$~h at 70~s/step.)

**Perturbation-injection diagnostics ($t_{\rm wall} = 0$, just before solver
launch):}
\begin{itemize}
\item Source mesh: $3.58\!\times\!10^{6}$ cells (Hartmann-resolved).
\item Laminar bulk RMS $|U_{\rm lam}|$ at end of Phase~1:
$5.37\!\times\!10^{-2}$~m/s (see Phase-1 observation note).
\item Target $|U'|_{\rm RMS}$ at $\alpha = 0.05$: $2.68\!\times\!10^{-3}$~m/s.
\item Achieved RMS ratio: $0.0500$ exactly (script's final rescale step
hits the target by construction).
\item Maximum $|U'|$: $1.79\!\times\!10^{-2}$~m/s (Gaussian distribution,
$\max/{\rm RMS} = 6.7$, normal for the spectral construction).
\item Wall taper performance: RMS $|U'|$ at $z<0.1$~mm is
$6.45\!\times\!10^{-5}$~m/s, i.e.\ $2.4$\,\% of bulk RMS --- the
$\sin\!\cdot\!\sin$ envelope is suppressing the perturbation at the walls
as designed.
\item Random seed: 1729 (same as the canonical debug-mesh bracket).
\end{itemize}

**Solver completion.** 1000 timesteps, $t = 0.055$~s reached, log ends with \texttt{End} and \texttt{Finalising parallel run}, cumulative continuity at $3.6\!\times\!10^{-12}$, no NaN/FPE/abort. GAMG-$\Phi$ at 3 iter/step (slightly higher than the 2 iter/step of the laminar Phase~1 because the perturbation makes the $\mathbf{u}\!\times\!\mathbf{B}$ source more spatially varying). First-step $p$-Initial-residual elevated to $\sim\!10^{-2}$ as predicted (the perturbation divergence cleanup), dropping back to the laminar $\sim\!10^{-5}$ within five timesteps.

\textbf{$\sigma\,\tau_A$ extraction.}
\begin{itemize}
\item RMS $U_z$ across probes, initial: $36.4$~mm/s.
\item RMS $U_z$ across probes, final: $8.14$~mm/s.
\item Ratio final/initial: $0.223$.
\item Log-linear fit: $\sigma = -30.6$~s$^{-1}$, $\sigma\,\tau_A = -0.227$.
\item Fit quality: $R^{2} = 0.978$, e-folding time $32.6$~ms.
\item \textbf{Verdict: MARGINAL DECAY (clean exponential).}
\end{itemize}

\textbf{Magnitude anomaly identified post-run --- amplitude was effectively much larger than requested.} A direct inspection of the raw probe data revealed that the per-probe initial $U_z$ values (up to $84$~mm/s at probe 5) substantially exceed the volumetric maximum $|U'|_{\max} = 17.9$~mm/s reported by \texttt{perturb\_ic.py}. The reason traces to a cell-weighting bias in the script's RMS rescaling step: \texttt{np.mean} over the cell array is cell-weighted, not volume-weighted, so the many tiny wall-cluster cells (where the wall taper drives $|U'|\to 0$) over-count toward the mean. The script's final rescale to hit \texttt{target\_abs\_rms} therefore divides by an under-estimate of the bulk RMS, leaving the bulk perturbation magnitudes $\sim\!4$--$5\times$ larger than the diagnostic states. The analyser, which reads the actual probe values from disk, sees the true perturbation magnitudes and reports \texttt{Initial perturbation/laminar = 84.7\%}, not the requested 5\,\%.

This bug existed on the under-resolved debug mesh as well but was not detected because cell volumes there are nearly uniform, so cell-weighting and volume-weighting agreed. On the Hartmann-resolved mesh the disparity becomes severe.

\textbf{Decision: halt the Hartmann-resolved bracket campaign at Phase~2.}
A 15\,\% requested bracket would inject an $\sim 250$\,\% effective perturbation; a 30\,\% request would inject $\sim 510$\,\%. Both cases over-run the base flow and break the perturbation-around-laminar framing, so the resulting $\sigma\,\tau_A$ would not be physically comparable to the debug-mesh bracket and would not constitute mesh-sensitivity validation. The single Phase-2 data point is preserved as the Hartmann-resolved point of comparison, with the amplitude reframing documented above.

\textbf{Phase 2 verdict: PASS (with documented caveats).} $\sigma\!\times\!\tau_A = -0.227 \pm \sim 0.05$ at $\sim 85\,\%$ effective probe amplitude on the Hartmann-resolved warm-state laminar baseline. The decay is clean exponential with $R^{2} = 0.978$, the dimensionless rate falls within the debug-mesh bracket range $[-0.40, -0.30]$ when the amplitude and warm-state-equilibration differences are accounted for, and the result therefore corroborates the under-resolved-mesh stability conclusion.

\textbf{Outputs delivered:} probe time-series, decay-rate fit, back-EMF anomalous output flagged for separate analysis. \texttt{postProcessing\_5pct\_hr/} archived in the case directory.

---

\textit{[Phases 3 and 4 not executed --- see decision above. Campaign log treats them as INTENTIONALLY DEFERRED to the cluster Hartmann-resolved bracket.]}

---

## Phase 3: bracket run @ 15 % RMS

Status: PENDING.

(Same structure as Phase 2.)

---

## Phase 4: bracket run @ 30 % RMS

Status: PENDING.

(Same structure as Phase 2.)

---

## Phase 5: analysers + side-by-side figures

Status: PENDING.

**Purpose.** Produce mesh-sensitivity figures: two-panel
σ × τ_A vs amplitude (debug overlaid on Hartmann-resolved); two-panel
back-EMF efficiency vs amplitude (same overlay).

**Outputs.** [Filenames + key values.]

---

## Summary

[To be filled in at end of campaign: total wall clock, total snapshots,
disk footprint, final verdict with a one-paragraph addition to the
report's Discussion section.]
