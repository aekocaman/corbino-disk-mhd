# corbino_turbulance_debug_hartmann — Hartmann-resolved laptop sibling case

A sibling of `corbino_turbulance_debug/`: same topology, same physics, same
solver, same protocol — **only** the z-grading is finer so that the Hartmann
boundary layer δ_H ≈ 53 μm is resolved with ≥ 5 cells (per the criterion in
finding #24).  Intended for an overnight laminar baseline + a tomorrow-after
perturbation bracket, to compare against the under-resolved `corbino_turbulance_debug/`
results and provide a mesh-sensitivity data point for the report.

> **Still not the production DNS mesh.**  r and θ resolutions are unchanged
> from the under-resolved debug case.  This mesh adds Hartmann-layer
> resolution only, which is the dimension most directly tied to the
> absolute U_∞ amplitude and the wall shear / torque closure.

## Sizing

  - Cells per block: 200 × 96 × 70  (vs debug 200 × 96 × 30, vs production 600 × 384 × 70).
  - Total: 4 × 200 × 96 × 70 = **5.376 × 10⁶ cells** (~2.3× larger than debug).
  - Working set: ~8 GB at ~1.5 KB/cell — borderline on a 16 GB laptop, fine if other apps are closed.
  - Decomposition: 8 Scotch subdomains → ~0.67 M cells / rank.

## z-grading (Hartmann-resolved)

  - Wall layer: `Δz_min = 7 μm`, `q = 1.2`, 11 cells/side.
  - Side span = 7 μm × (1.2¹¹ − 1) / 0.2 ≈ 225 μm.
  - Cells inside δ_H = 53 μm: **8** (cells 1–8 cumulative widths 7…52 μm).
    "≥ 5 cells in δ_H" criterion satisfied.
  - Bulk: 48 uniform cells × 199 μm.
  - Total: 11 + 48 + 11 = 70.

## Compute budget (vs under-resolved debug)

| metric                       | debug (1 cell in δ_H) | hartmann (8 cells) |
|------------------------------|-----------------------|--------------------|
| total cells                  | 2.30 M                | 5.38 M             |
| Δt cap (CFL@dz_min, U=0.1)   | 5e-4 s                | 3.5e-5 s           |
| 0 → 0.5 s laminar baseline   | ~50 min               | ~7-9 h             |
| 50 ms perturbation restart   | ~3 min                | ~50 min            |
| 50 snapshot disk             | ~4.5 GB binary        | ~10.5 GB binary    |

## Setup procedure (WSL, mirrors the debug case)

```bash
# from ~/OpenFOAM/reybo78-13/run/ :
cp -r corbino_turbulance_debug corbino_turbulance_debug_hartmann
cd corbino_turbulance_debug_hartmann

# Stage the Hartmann-resolved mesh + dicts
PROJ="/mnt/c/Users/user/Documents/Claude/Projects/Turbulance Instabilities of a conducting fluid in 3-D Corbino Geomtry"
cp "$PROJ/corbino_turbulance_debug_hartmann/system/blockMeshDict"     system/blockMeshDict
cp "$PROJ/corbino_turbulance_debug_hartmann/system/decomposeParDict"  system/decomposeParDict
cp "$PROJ/corbino_turbulance_debug_hartmann/system/controlDict.long"  system/controlDict.long
# Re-use the GAMG-Phi fvSolution from the debug case (already verified):
cp "$PROJ/corbino_turbulance_debug/system/fvSolution.long"            system/fvSolution.long

# Purge inherited mesh, run state, perturbation backups, and outputs
rm -rf 0/  0.* constant/polyMesh/  processor*  postProcessing*  log.* epotmhd.pid
cp -r ../corbino_turbulance_debug/0  0   # restore IC templates (will be remapped)

# Build the new mesh
blockMesh                                  # ~30-60 s, 5.38 M cells
checkMesh -allTopology -allGeometry > log.checkMesh
tail -30 log.checkMesh                     # confirm "Mesh OK"
```

`checkMesh` should report:

  - 5,376,000 cells, 4 hex blocks
  - max non-orthogonality < 60° (target)
  - max skewness < 4 (target)
  - max aspect ratio ~ 240  (driven by `r_out · Δθ / Δz_min` — bigger than debug
    because Δz_min is 7× smaller, but still inside OpenFOAM's tolerance)

## Overnight launch (the laminar baseline)

```bash
cd ~/OpenFOAM/reybo78-13/run/corbino_turbulance_debug_hartmann

# Stage the long-run dicts
cp system/controlDict.long  system/controlDict
cp system/fvSolution.long   system/fvSolution

# Decompose
decomposePar > log.decomposePar.long 2>&1
tail -3 log.decomposePar.long

# Launch detached so it survives the laptop sleeping (caveat: not all
# laptops keep WSL processes alive across full sleep; check 'powercfg
# /requests' on the Windows side or run with the lid open)
nohup mpirun -np 8 epotMhdFoam -parallel > log.epotMhdFoam.long 2>&1 &
echo $! > epotmhd.pid
disown
echo "Launched PID $(cat epotmhd.pid).  Expect ~7-9 h wall clock."
```

## Morning health check

```bash
cd ~/OpenFOAM/reybo78-13/run/corbino_turbulance_debug_hartmann

# Still alive?  If 'ps' returns the row, it's running; if no row, exited.
ps -p $(cat epotmhd.pid) -o pid,etime,pcpu,pmem,cmd

# Progress
grep -c "^Time = " log.epotMhdFoam.long
tail -40 log.epotMhdFoam.long
ls processor0/ | grep -E '^[0-9]' | sort -g | tail -10

# If 'End' is the last log line, done.  If still running, check its rate:
#   ETA  =  (0.5 - latest_time) / (latest_time / elapsed_seconds)
```

## See also

  - `../corbino_turbulance_debug/README.md` — the under-resolved sibling
  - `../report/step1_laminar_report.md` — finding #24 (mesh design),
    #26 (under-resolved-mesh price), #28-30 (Step-2 first results)
