# Step 1 — 360° Laminar Corbino Baseline

Full-annulus counterpart to the 90° wedge case from the progress report.
Run this before turbulence (Step 2); the time-averaged azimuthal profile
must match the wedge result within the performance targets from section
*Model Parameters and Performance Goals*.

---

## What changed vs. the wedge case

| Area | Wedge (old) | 360° (this case) | Why |
|---|---|---|---|
| Geometry | 90° sector with `sides` symmetry patches | Full annulus, 4 hex blocks | Turbulent fluctuations (Step 2) break the azimuthal symmetry; the wedge is no longer a valid simplification |
| Patches | `innerElectrode`, `outerElectrode`, `insulatingPlates`, `sides` (symmetry) | `innerElectrode`, `outerElectrode`, `insulatingPlates` | No symmetry plane in 360°; removing `sides` is what the BC files encode |
| `deltaT` | Fixed 5 × 10⁻³ s | Adaptive, `maxCo 0.5`, starts at 10⁻⁷ s | Wedge run hit Courant = 26,880 at t = 0.035 s → sigFpe. See `controlDict` header |
| Inner `Phi` | `uniform 1.0` (≈ 135 kA!) | Linearly ramped 0 → 3.703 × 10⁻³ V over 50 ms | Physical V for I = 500 A; ramp over ~7·τ_A avoids the impulsive Lorentz kick |
| Mesh | 153,600 cells, 90° wedge | 614,400 cells (60 × 80 × 32 per block × 4) | Same per-radian resolution as wedge; reduce `nθ` to 40 for a faster first run |

---

## Directory layout

```
case_360_laminar/
├── 0/
│   ├── B              # B = 0.5 T zhat, fixedValue on all walls
│   ├── Phi            # ramped electrode potential
│   ├── U              # no-slip everywhere
│   ├── p              # fixedFluxPressure; closed domain, see note below
│   └── pB             # magnetic-pressure projection, zeroGradient
├── constant/
│   └── transportProperties   # Galinstan σ, ν, ρ, µ
├── system/
│   ├── blockMeshDict         # 4-block O-grid, hyperbolic z-grading
│   ├── controlDict           # adaptive Δt, function objects
│   ├── fvSchemes             # 2nd-order linear, linearUpwind convection
│   └── fvSolution            # GAMG for p, PCG for Φ, pBiCGStab for U/B
└── README.md          (this file)
```

---

## Run sequence (WSL)

```bash
cd "case_360_laminar"

# 1. Generate the mesh
blockMesh 2>&1 | tee log.blockMesh

# 2. Sanity-check topology
checkMesh 2>&1 | tee log.checkMesh
# Expect:  Mesh OK.
# Watch for:
#   - Max aspect ratio (target < 1000; ~800 is fine for our grading)
#   - Max non-orthogonality (keep < 65°; bump nNonOrthogonalCorrectors if ~70°)
#   - Negative-volume cells => vertex-order bug, see "Small Determinant"
#     appendix in the progress report

# 3. Run the solver
mhdFoam 2>&1 | tee log.mhdFoam
```

If `checkMesh` complains about **severe non-orthogonality** at the
block-to-block interfaces, it's almost always the arc edges disagreeing on
the midpoint side.  All eight arcs above use the **+45° convention**: each
arc mid-point is at (±R/√2, ±R/√2), chosen so the arc bulges *outward*
(larger radius) in every quadrant.  If you reverse even one, you get the
same negative Jacobian you hit on the wedge.

---

## Validation targets (same as progress report)

| Metric | Definition | Target |
|---|---|---|
| `E_τ` | \|τ_num − τ_ana\| / τ_ana,  τ_ana = I B (r_out² − r_in²) / 2 | ≤ 0.01 |
| HLR | δ_Ha / Δz_wall,  δ_Ha ≈ h / Ha | ≥ 8 |
| RMSE_Φ | RMS of Φ − Φ_theory where Φ_theory ∝ 1/r in core | ≤ 0.015 |
| Rm | µ₀ σ U₀ L | < 0.1 |
| ε_sym | \|f_r\| / \|f_θ\| | < 10⁻³ |

**Pass criterion for this step:** all five targets met *and* the
time-averaged u_θ(r, z) profile on any θ = const slice matches the wedge
result to within the mesh-convergence bound (|ΔU| / U_0 < 0.5 %).

---

## What to change first if it blows up

1. **CFL still runaway?** Lengthen the ramp in `0/Phi` to 0.1 s (table value 0.1
   instead of 0.05).  Each doubling of ramp time halves the peak du_θ/dt.
2. **Hartmann layer under-resolved after `checkMesh`?**  Drop `nz = 32` to
   `nz = 40` and adjust grading so the first cell stays ≤ δ_Ha / 8.
3. **pBiCGStab stalling on U?**  Tighten `nNonOrthogonalCorrectors` to 2.
4. **Torque error > 1 %?**  First suspect is the mesh — run a Richardson
   extrapolation with `(40 60 24)` and `(60 80 32)` before blaming the BCs.

---

## Post-processing skeleton

Once a snapshot is written:

```bash
postProcess -func "writeCellCentres"    # exports cell coordinates for RMSE_Φ
postProcess -func "mag(U)"              # speed field for Courant audit
postProcess -func "components(U)"       # Ux, Uy, Uz so you can build u_θ

# τ_num via the cross-product integrand (progress report Eq. 17):
#    ∫ (x * f_y − y * f_x) dV
# Compute f = J × B in ParaView's Calculator then Volume Integrate.
```

τ_ana with the design parameters:

```
τ_ana = I B (r_out² − r_in²) / 2
      = 500 × 0.5 × (0.1² − 0.02²) / 2
      = 500 × 0.5 × 0.0096 / 2
      = 1.20 N·m
```

So a passing run reports τ_num in the range **1.188 … 1.212 N·m**.
