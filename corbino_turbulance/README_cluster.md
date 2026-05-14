# Cluster runbook — corbino_turbulance Step-2 production DNS

Start-to-finish instructions for queuing the 64.5 M-cell DNS on a
TÜBİTAK ULAKBİM TRUBA / Bilkent CSE cluster. Self-contained: anyone
with cluster credentials and basic SLURM familiarity should be able to
follow this top-to-bottom and produce the headline result.

**Authoritative source for every parameter and rationale:**
`report/step1_laminar_report.md` (the consolidated technical report).
This README is a *deployment* document, not a methodology document.

---

## 0. Preflight — read before submitting anything

  1. Confirm cluster supports OpenFOAM 13 Foundation release. **Not** ESI
     v2312 / v2406 — our custom solver `epotMhdFoam` is binary-incompatible
     across Foundation/ESI variants.
     ```
     module avail openfoam        # check what's on the cluster
     ```
     If only ESI is available, file a ticket with the cluster admin
     to install OpenFOAM 13 Foundation, or build it in your home
     directory (~30 min on a dedicated build node, instructions in
     `epotMhdFoam/README.md`).
  2. Confirm aggregate RAM availability. Working set is ~120 GB at
     256 ranks (~470 MB/rank). Allow 50 % headroom →
     **≥ 192 GB aggregate** must be reachable.
  3. Confirm scratch-disk quota. One full Phase-5 run produces:
     - 90 snapshots × 2 GB binary    ≈ 180 GB
     - probe data                     ≈ 1 GB
     - logs and decomposition maps    ≈ 5 GB
     **Allocate ≥ 250 GB on scratch** before submission.
  4. Confirm walltime allowance. One Phase-5 run at 256 ranks is
     ≈ 50–80 hours of wall clock. Most TRUBA partitions cap individual
     jobs at 72 h; if the cluster's `long` partition allows 168 h, use
     that. Otherwise plan for **two consecutive 72 h jobs** with
     `startFrom latestTime;` chaining.
  5. Verify that you have the email address right in `submit_corbino.sh`
     (`--mail-user=...`). The job will email you at BEGIN, END, FAIL.

---

## 1. Stage the case onto cluster scratch

```bash
# 1a. Get the case files from your laptop (or git remote)
cd $SCRATCH                                     # cluster scratch directory
mkdir -p corbino_dns && cd corbino_dns

# Option A: rsync from a logged-in laptop (most common)
rsync -avz --exclude 'processor*' --exclude '*.tar.gz' \
    eraylaptop:/path/to/corbino_turbulance/  ./

# Option B: scp the tar archive that travels with the report bundle
scp eraylaptop:/path/to/corbino_turbulance.tar.gz .
tar -xzf corbino_turbulance.tar.gz && rm corbino_turbulance.tar.gz

# 1b. Build the custom solver (one-shot per cluster, ~10 min)
mkdir -p $HOME/OpenFOAM/$USER-13/applications/solvers
cd       $HOME/OpenFOAM/$USER-13/applications/solvers
rsync -avz $SCRATCH/corbino_dns/../epotMhdFoam ./       # source tree
cd epotMhdFoam
source /truba/sw/centos7.9/comp/openfoam/13/etc/bashrc  # adjust path
wmake
which epotMhdFoam                  # confirm: $FOAM_USER_APPBIN/epotMhdFoam
```

If `wmake` fails: see `epotMhdFoam/README.md` for the three known
build pitfalls (missing header path, missing -lfoam library, fragile
`fvm::Sp` overload).

---

## 2. Initial-condition warm-start via mapFields

There are two ways to give the production solver a sensible IC:

  **(a) Cold start from rest.** Set `0/U` to zero everywhere except
       boundary patches; ramp ΔV from 0 → 3.7 mV over 50 ms (already
       configured in `0/Phi`). Cost: ≈ 1 τ_H of bulk equilibration plus
       ≈ 2 τ_H of Hartmann-layer relaxation = 3 τ_H ≈ 2 s of wasted
       physical time before averaging window opens.

  **(b) Warm start via mapFields from the laptop debug case** (recommended).
       Copy the converged laminar field at $t=0.5$ s from
       `corbino_turbulance_debug/processor0/0.5/` (single-domain
       reconstructed) and `mapFields` it onto the production mesh.
       Saves ≈ 1.5 s of wasted compute per run.
       ```bash
       # On the cluster, in $SCRATCH/corbino_dns/:
       DEBUG_CASE=$SCRATCH/corbino_dns_debug    # path to the staged debug case
       mapFields "$DEBUG_CASE" -sourceTime 0.5 -consistent \
           > log.mapFields 2>&1
       tail -5 log.mapFields
       ls 0/                                    # confirm U, Phi, p, B updated
       ```

The 17 % U_∞ shortfall in the debug-mesh laminar (finding #26 of the
report) is inherited; the warm-started run still needs ≥ 1 τ_H of bulk
equilibration to recover the true Hartmann-resolved $U_{\infty}$, but
this is much cheaper than starting from rest.

---

## 3. Submit the job

```bash
cd $SCRATCH/corbino_dns

# Edit submit_corbino.sh:
#   - --account, --partition, --qos to match your cluster
#   - --mail-user to your real email
#   - --ntasks (default 256; can scale 128/512 — must match decomposeParDict)
#   - --time matches partition cap (72 h typical, 168 h on long QoS)

vim submit_corbino.sh        # or your favourite editor

# Submit:
sbatch submit_corbino.sh

# Track:
squeue -u $USER                                # queue position
sacct -j <job_id> --format=jobid,state,elapsed # post-completion summary
```

---

## 4. Monitor a running job

```bash
# Tail the live log (refresh ~every 10 s):
tail -f $SCRATCH/corbino_dns/log.epotMhdFoam.<JOB_ID>

# Solver progress:
grep -c "^Time = " log.epotMhdFoam.<JOB_ID>
grep "ExecutionTime" log.epotMhdFoam.<JOB_ID> | tail -3
grep "cumulative" log.epotMhdFoam.<JOB_ID> | tail -3
grep "GAMG.*Phi," log.epotMhdFoam.<JOB_ID> | awk '{print $NF}' | sort -n | tail -5

# Snapshots written so far:
ls processor0/ | grep -E '^[0-9]' | sort -g | tail -10
```

**Healthy signs:**
  - `Phi GAMG` 2–10 iterations/step
  - `cumulative` continuity error < 1e-9 by the end
  - `ExecutionTime` increasing steadily (~1 s/step at 256 ranks; if much
    larger, check for swap pressure or interconnect issues)
  - No `FATAL`, `nan`, `inf`, `exception` strings in the log

**If something goes wrong:**
  - Job killed for time-limit overrun → resubmit with
    `startFrom latestTime;` already in `controlDict.production`,
    no other change needed.
  - Solver aborts with `nan` → almost always indicates Δt collapse.
    Edit `system/controlDict` while the next job is queued: tighten
    `maxCo` from 0.3 to 0.2, restart from the last good snapshot.
  - Disk full → snapshots are still on `processor*/`; reduce
    `writeInterval` (e.g. 0.1 s instead of 0.05 s) and resume.

---

## 5. Post-processing on the cluster

`reconstructPar` is required to produce single-domain snapshots that
the analyser scripts expect. It is expensive at 64.5 M cells (~30 min
per snapshot, single-threaded) and cluster-scratch-bound.

```bash
# Reconstruct only the snapshots you need for analysis (e.g. last 50):
LAST_SNAPSHOTS=$(ls processor0/ | grep -E '^[0-9]' | sort -g | tail -50)
for t in $LAST_SNAPSHOTS; do
    reconstructPar -time $t > log.reconstructPar.$t 2>&1
done

# Then run the analyser scripts (laptop scripts, work on the cluster):
foamFormatConvert -time <t>           # binary → ASCII for Python parsers
postProcess -func writeCellCentres -time <t>
python3 ../scripts/perturb_ic.py <t> --rms 0.05         # if perturbation run
python3 ../scripts/analyse_perturb.py --probes-dir postProcessing/probes/<t>
python3 ../scripts/analyse_backemf.py  --probes-dir postProcessing/probes/<t>
python3 ../scripts/make_bracket_report.py
```

**Recommended workflow:** keep `processor*/` directories on cluster
scratch; copy only `postProcessing/`, the reconstructed snapshots,
and `report/` figures back to your laptop for analysis. This saves
~150 GB of wide-area transfer per run.

---

## 6. Production sweep schedule (staged, 4–6 weeks total)

Per the eight-phase plan in `report/step1_laminar_report.md` §9.2:

| Phase | Description | Cluster time |
|-------|-------------|--------------|
| Phase 2 | Refined-mesh laminar verification (1 τ_H) | ~50 h on 64 ranks |
| Phase 3 | Perturbed IC injection + smoke test | ~3 h on 256 ranks |
| Phase 4 | Spin-up + go/no-go (1 τ_H) | ~50 h on 256 ranks |
| Phase 5 | Statistically-stationary 5 τ_H run | ~150 h on 256 ranks |
| Phase 6 | Diagnostics (post-processing only) | local + cluster, ~1 week |
| Phase 7 | GCI: coarse / medium / fine | ~250 h on 256 ranks (over 3 meshes) |
| Phase 8 | Manuscript + Zenodo release | local |

**Total cluster compute estimated 500–600 wall-hours at 256 ranks ≈
130 000 core-hours.** Allocate accordingly in your TRUBA proposal.

---

## See also

  - `preflight_checklist.md` — single-page sanity sheet for the cluster admin
  - `submit_corbino.sh` — the SLURM template referenced from §3
  - `system/{controlDict,fvSolution,fvSchemes,decomposeParDict}.production`
    — the production-tuned case files
  - `report/step1_laminar_report.md` — methodology and prior results
  - `report/audit_2026_05_03.md` — sanity / unit / literature audit
