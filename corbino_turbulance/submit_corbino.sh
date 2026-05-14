#!/bin/bash
# =====================================================================
#  SLURM submit script for the Corbino DNS production run
#
#  Targets TÜBİTAK ULAKBİM TRUBA (and Bilkent CSE) cluster partitions.
#  Adjust --partition, --account, --qos and module-load lines per the
#  cluster's local conventions before first submit.
#
#  Usage:
#      sbatch submit_corbino.sh
#  or with overrides:
#      sbatch --ntasks=512 --time=72:00:00 submit_corbino.sh
#
#  Pre-flight: see preflight_checklist.md for required modules,
#  expected RAM/disk, and sanity checks the cluster admin should run.
# =====================================================================

#SBATCH --job-name=corbino-dns
#SBATCH --account=YOUR_ACCOUNT          # <-- replace with your TRUBA / CSE account code
#SBATCH --partition=long                # <-- TRUBA: "long" or "compute"; CSE: "common"
#SBATCH --qos=long                      # <-- match partition QoS

#SBATCH --ntasks=256                    # MPI ranks; matches decomposeParDict
#SBATCH --nodes=8                       # 32 ranks/node typical (Intel Xeon Gold)
#SBATCH --ntasks-per-node=32

#SBATCH --time=72:00:00                 # 72 h = 3 days; one full Phase-5 production run
#SBATCH --mem=512G                      # 64 GB/node x 8 nodes; 8 GB/rank average

#SBATCH --output=slurm-%j.out           # %j = job ID
#SBATCH --error=slurm-%j.err

#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ahmeteraykocaman@gmail.com   # <-- check, then change if needed

# ---------------------------------------------------------------------
# Environment: load OpenFOAM 13 and OpenMPI
# ---------------------------------------------------------------------
# TRUBA-style modules (verify on your cluster):
module purge
module load openmpi/4.1.4-gcc-9.3.0
module load openfoam/13         # Foundation release; if cluster only has
                                # ESI v2312, contact admin -- our solver
                                # epotMhdFoam was built against Foundation
                                # 13 and is binary-incompatible with ESI

# Source OpenFOAM 13 environment (adjust path to your cluster's install):
source /truba/sw/centos7.9/comp/openfoam/13/etc/bashrc

# Sanity prints to slurm-%j.out
echo "===== Run metadata ====="
echo "Job ID:      $SLURM_JOB_ID"
echo "Job name:    $SLURM_JOB_NAME"
echo "Nodes:       $SLURM_JOB_NUM_NODES"
echo "Tasks:       $SLURM_NTASKS"
echo "Submit dir:  $SLURM_SUBMIT_DIR"
echo "Date:        $(date -Iseconds)"
echo "OpenFOAM:    $WM_PROJECT_DIR ($(which epotMhdFoam))"
echo

cd "$SLURM_SUBMIT_DIR" || exit 1

# ---------------------------------------------------------------------
# Stage production case files (overwrite the default names)
# ---------------------------------------------------------------------
# These three are the production-tuned dictionaries banked in repo:
cp system/controlDict.production  system/controlDict
cp system/fvSolution.production   system/fvSolution
cp system/fvSchemes.production    system/fvSchemes

# decomposeParDict has numberOfSubdomains hard-coded to 64 by default;
# scale to match $SLURM_NTASKS:
sed -i "s/^numberOfSubdomains.*/numberOfSubdomains $SLURM_NTASKS;/" \
    system/decomposeParDict
echo "decomposeParDict set to $SLURM_NTASKS subdomains."

# ---------------------------------------------------------------------
# Mesh + decomposition (cheap; rebuild every job to be safe)
# ---------------------------------------------------------------------
# Skip if mesh already exists and is current:
if [ ! -d constant/polyMesh ]; then
    echo "===== blockMesh ====="
    blockMesh > log.blockMesh 2>&1
    tail -3 log.blockMesh
fi

if [ -d processor0 ] && \
   [ "$(ls processor0/ 2>/dev/null | head -1)" = "constant" ]; then
    echo "Decomposition exists; skipping decomposePar."
else
    echo "===== decomposePar ====="
    decomposePar -force > log.decomposePar 2>&1
    tail -3 log.decomposePar
fi

# ---------------------------------------------------------------------
# Run the solver
# ---------------------------------------------------------------------
echo "===== Launching epotMhdFoam on $SLURM_NTASKS ranks ====="
echo "Start: $(date -Iseconds)"

# Use srun (preferred on SLURM) rather than mpirun; respects cgroup
# limits and hands off correctly to the queue manager:
srun --mpi=pmi2 epotMhdFoam -parallel \
     > log.epotMhdFoam.${SLURM_JOB_ID} 2>&1

EXIT=$?
echo "End:   $(date -Iseconds)"
echo "Solver exit code: $EXIT"

# ---------------------------------------------------------------------
# Post-run summary (cheap; tails the log even on failure)
# ---------------------------------------------------------------------
echo "===== Post-run summary ====="
echo "Steps completed: $(grep -c '^Time = ' log.epotMhdFoam.${SLURM_JOB_ID})"
echo "Final continuity errors:"
grep "cumulative" log.epotMhdFoam.${SLURM_JOB_ID} | tail -3
echo
echo "Final ExecutionTime:"
grep "ExecutionTime" log.epotMhdFoam.${SLURM_JOB_ID} | tail -1
echo

if [ $EXIT -eq 0 ]; then
    echo "Run finished successfully."
    # Optional: kick off reconstructPar on every snapshot (heavy; ~30 min/snapshot
    # at 64.5 M cells; consider a separate dependent job).
    # sbatch --dependency=afterok:$SLURM_JOB_ID reconstruct_corbino.sh
else
    echo "Run aborted; do not delete processor*/ directories."
fi

exit $EXIT
