# corbino_turbulance_debug — laptop wiring case

This is a **deliberately under-resolved** sibling of `corbino_turbulance/`,
intended to run on a laptop (≤16 GB RAM) for the sole purpose of
shaking down the case wiring before the production DNS is queued on a
cluster.

> **Not a physics-grade mesh.**  ~1 cell across `δ_H ≈ 53 µm`, vs the DNS
> requirement of ≥5 cells.  

## It is *for*

  - Verifying `blockMesh` / `checkMesh` topology and patch identity.
  - Exercising `decomposePar` (Scotch, 8 ranks).
  - Smoke-testing `mhdFoam` (PIMPLE start-up, ddt/divSchemes choices).
  - Validating the divergence-free turbulent IC injection on `0/U`
    (Phase 3 of `report/step1_laminar_report.md` §9).
  - Shaking down snapshot I/O cadence and `runTimeModifiable` controls.
  - Smoke-testing the post-processing pipeline (paraview / pyFoam /
    FFT scripts).

## It is *not* for

  - Quantitative spectra `E(k_θ)`, Reynolds stresses, dissipation budgets.
  - GCI / mesh-convergence study (use the coarse/medium/fine variants
    in `corbino_turbulance/`).
  - Any publishable  result 

## Sizing

  - Cells per block: 200 × 96 × 30  (vs production 600 × 384 × 70).
  - Total: 4 × 200 × 96 × 30 = **2.304 × 10⁶ cells** (~28× smaller).


## z-grading (relaxed)

  - Wall layer: `Δz_min = 50 µm`, `q = 1.3`, 5 cells/side.
  - Bulk: 20 uniform cells × ~455 µm.
  - Total: 5 + 20 + 5 = 30, side span = 452 µm, bulk span = 9.096 mm.

## r- and θ-grading (matched to production)

  - r: 1.05 : 1 from inner to outer electrode.
  - θ: 1.0 (uniform; FFT-friendly).

## Provenance

Topology, vertices, edges, and patches are inherited verbatim from
`../corbino_turbulance/system/blockMeshDict`.  Only the cell counts and
the z-grading multi-segment list differ.

## Setup procedure (WSL, mirrors production)

```
# from ~/OpenFOAM/reybo78-13/run/ :
cp -r case_360_laminar corbino_turbulance_debug
cd corbino_turbulance_debug

# overwrite the mesh and decomposition dicts with the debug versions:
cp /mnt/c/Users/user/Documents/Claude/Projects/'Turbulance Instabilities of a conducting fluid in 3-D Corbino Geomtry'/corbino_turbulance_debug/system/blockMeshDict   system/blockMeshDict
cp /mnt/c/Users/user/Documents/Claude/Projects/'Turbulance Instabilities of a conducting fluid in 3-D Corbino Geomtry'/corbino_turbulance_debug/system/decomposeParDict system/decomposeParDict

rm -rf 0/ 0.* constant/polyMesh/   # purge laminar fields & old mesh
cp -r ../case_360_laminar/0 0      # restore IC templates (will be remapped)
blockMesh                          # rebuild mesh (~10-20 s, 2.3 M cells)
checkMesh -allTopology -allGeometry > log.checkMesh
```

`checkMesh` should report:

  - 2,304,000 cells, 4 hex blocks
  - max non-orthogonality < 60° (target)
  - max skewness < 4 (target)
  - max aspect ratio < ~35 (driven by `r_out · Δθ / Δz_min`)

## See also

`../corbino_turbulance/README.md` for the production-DNS sibling, and
`report/step1_laminar_report.md` §8–§9 for the turbulence roadmap and
eight-phase execution plan.
