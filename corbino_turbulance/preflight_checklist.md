# Pre-flight checklist — Corbino DNS production deployment

**One-page summary for the cluster administrator and PI.**
For the full deployment recipe see `README_cluster.md`.

## Software requirements

| component | version | notes |
|---|---|---|
| Operating system | Linux (CentOS 7+, Rocky 8+, Ubuntu 22.04) | any modern HPC distro |
| MPI | OpenMPI 4.x or MPICH 3.4+ | tested on OpenMPI 4.1.4 |
| OpenFOAM | **13 (Foundation release)** | `https://openfoam.org/release/13/`. **Not** ESI v2312/v2406 — solver is binary-incompatible across forks |
| Compiler | GCC 9.3+ | needed for OpenFOAM 13 build |
| Python (post-proc) | 3.10+ | NumPy, matplotlib; on login node only |

## Hardware footprint per production job

| resource | requirement | comment |
|---|---|---|
| MPI ranks | 128–512 | optimum 256; weak-scales near linearly to 512 |
| RAM/rank | ~470 MB at 256 ranks; ~940 MB at 128 | 8 GB/rank is comfortable |
| Aggregate RAM | ≥ 192 GB | mesh + buffers + 2×snapshot |
| Interconnect | InfiniBand or equivalent | pressure-Poisson all-reduces are latency-sensitive |
| Wall clock per Phase-5 run | 50–80 h at 256 ranks | one job at standard 72 h walltime, or two chained 36 h jobs |
| Scratch disk per run | ~250 GB | snapshots dominate (90 × 2 GB binary) |

## Check 1: solver builds cleanly on cluster

```bash
module purge
module load gcc/9.3.0 openmpi/4.1.4 openfoam/13
source $(spack location -i openfoam)/etc/bashrc

cd $HOME/OpenFOAM/$USER-13/applications/solvers/epotMhdFoam
wmake
which epotMhdFoam      # should resolve to $FOAM_USER_APPBIN/epotMhdFoam
```

Expected output: a binary at `$FOAM_USER_APPBIN/epotMhdFoam`, ~2 MB.
If `wmake` fails, see `epotMhdFoam/README.md` for the three known
build pitfalls (paths, `-lfoam` linker, `fvm::Sp` overload).

## Check 2: `blockMesh` and `checkMesh` pass on a login node

```bash
cd $SCRATCH/corbino_dns
blockMesh > log.blockMesh 2>&1
tail -3 log.blockMesh   # expect "End"
checkMesh -allTopology -allGeometry > log.checkMesh 2>&1
grep -E "(cells:|hexahedra|Max non-orth|Max skewness|Max aspect|Mesh OK)" \
    log.checkMesh
```

Expected:
  - 64,512,000 cells, 4 hex blocks
  - Max non-orthogonality < 60°
  - Max skewness < 4
  - Max aspect ratio < 100 (production mesh has tighter ratio than the
    laptop debug case because Nθ = 384 is larger)
  - "Mesh OK" *or* "Failed 1 mesh checks" with the small-determinant
    warning only — see audit `report/audit_2026_05_03.md` Part 2 for
    why this warning is benign on Hartmann-resolved meshes.

`blockMesh` takes ~1–2 minutes; `checkMesh` ~5 minutes. Run on a login
node before submitting a multi-thousand-core-hour SLURM job.

## Check 3: tiny smoke test on 8–16 ranks

Before committing a 50-hour job, verify that the case advances on the
cluster's MPI stack. Recommended: a 100-step shadow job on 16 ranks,
~30 minutes wall clock at 64.5 M cells (~6× slower than 256 ranks but
trivially short).

```bash
# Make a copy of the case with reduced endTime:
cp -r $SCRATCH/corbino_dns $SCRATCH/corbino_smoke
cd $SCRATCH/corbino_smoke
sed -i 's/^endTime.*/endTime         5e-3;/'    system/controlDict.production
sed -i 's/^numberOfSubdomains.*/numberOfSubdomains 16;/' system/decomposeParDict
cp system/controlDict.production system/controlDict
cp system/fvSolution.production  system/fvSolution

# Submit a 1-hour job:
sbatch --time=01:00:00 --ntasks=16 --nodes=1 \
       --mem=128G submit_corbino.sh
```

If this completes with `End` in the log and continuity errors at
machine precision, the production submission is greenlit.

## Check 4: storage and quota

```bash
quota -s        # check user disk quota
df -h $SCRATCH  # check scratch capacity
```

Allocate at least 1 TB of scratch per production sweep (Phase 5 + 7
combined). The `processor*/` directories must remain on scratch
throughout — do not `tar` them as the run progresses.

## Email & notifications

The `submit_corbino.sh` template emails BEGIN, END, FAIL events.
Verify the email address matches your active inbox before first
submission.

---

**Sign-off:** when all four checks pass, the case is green for the
production sweep documented in `README_cluster.md` §6.
