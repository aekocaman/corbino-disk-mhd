# corbino_turbulance — Step-2 turbulent DNS case

This case directory holds the Step-2 turbulent DNS at
Re_a / Ha² ≈ 1.2 (PK14 Q2D-3D transition), full 360°.
Geometry, electrodes, and physical properties are inherited
verbatim from `case_360_laminar/`. Only the mesh, snapshot
cadence, and (later) the initial condition will differ.

## Provenance

The intended setup procedure (run by the user in WSL):

```
# from ~/OpenFOAM/reybo78-13/run/ :
cp -r case_360_laminar corbino_turbulance
cd corbino_turbulance

# overwrite the mesh dict with the DNS-graded version drafted
# in the project workspace:
cp /mnt/c/Users/user/Documents/Claude/Projects/'Turbulance Instabilities of a conducting fluid in 3-D Corbino Geomtry'/corbino_turbulance/system/blockMeshDict system/blockMeshDict
cp /mnt/c/Users/user/Documents/Claude/Projects/'Turbulance Instabilities of a conducting fluid in 3-D Corbino Geomtry'/corbino_turbulance/system/decomposeParDict system/decomposeParDict

rm -rf 0/ 0.* constant/polyMesh/   # purge laminar fields & old mesh
cp -r ../case_360_laminar/0 0      # restore IC templates (will be remapped)
blockMesh                          # rebuild mesh (~1-2 min, 64.5 M cells)
checkMesh -allTopology -allGeometry > log.checkMesh
```

`checkMesh` should report:

  - 64,512,000 cells, 4 hex blocks
  - max non-orthogonality < 60° (target)
  - max skewness < 4 (target)
  - max aspect ratio z-resolved < 50 (target)
  - max aspect ratio in θ-only ≈ 58 (acceptable; see §9.2 work record)

## Files (Step-2 specific)

  - `system/blockMeshDict` — (600, 384, 70) per block; 3-segment
    z-grading with Δz_min = 7 μm at each Hartmann wall (5 cells
    inside δ_H); 1.05:1 radial grading; uniform azimuthal.
  - `system/decomposeParDict` — Scotch decomposition, 64 → 128/256
    subdomains for production.
  - (TBD) `system/controlDict` — adaptive Δt with maxCo = 0.4,
    snapshot cadence 50 ms (Phase 5), total integration 6 τ_H ≈ 4.26 s.
  - (TBD) `0.org/U` — laminar-mapped + divergence-free turbulent
    perturbation, A_k_θ ∝ k_θ^(-7/6) (Phase 3).

## See also

`report/step1_laminar_report.md` §8 (turbulence roadmap) and
§9 (eight-phase execution plan, with per-phase math/values).
