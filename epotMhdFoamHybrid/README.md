# epotMhdFoamHybrid

Low-to-moderate magnetic-Reynolds-number incompressible MHD solver for OpenFOAM 13.

## Overview

`epotMhdFoamHybrid` combines two existing solvers:

- **`epotMhdFoam`** — the low-Rm electric-potential formulation, in which the induced magnetic field is neglected and current density is computed from the imposed `Phi` field via Ohm's law `J = sigma (-grad(Phi) + U x B)`.
- **Stock `mhdFoam`** — the magnetic-induction equation, which evolves the magnetic field `B` itself together with a solenoidal projection step.

This makes the solver usable across the full range from low Rm (where the induction equation produces essentially nothing) up to moderate Rm ~ O(1), where the induced `B'` becomes a finite physical contribution.

## Equations solved (per timestep)

The time loop applies operator splitting in the following order:

1. **Electric-potential Poisson** (uses current total `B`):

   ```
   div(sigma grad(Phi)) = div(sigma U x B)
   ```

2. **Lorentz body force** (decomposed for stability):

   ```
   F_L / rho = (sigma / rho) [ -grad(Phi) x B + (U.B) B - |B|^2 U ]
   ```

   The `-|B|^2 U` Hartmann-damping term is integrated implicitly via `fvm::Sp`, which is essential at high Hartmann number.

3. **Momentum predictor + PISO** for `U` and `p`.

4. **Magnetic induction equation** for the total field `B`:

   ```
   dB/dt + div(phi, B) - div(phiB, U) - laplacian(D_B, B) = 0
   ```

   with magnetic diffusivity `D_B = 1 / (mu sigma)`.

5. **Solenoidal projection** via magnetic-pressure scalar `pB`:

   ```
   laplacian(pB) = div(phiB)
   phiB := phiB - flux(grad(pB))
   B    := B    - grad(pB)
   ```

   This keeps `div(B)` bounded at the projection-tolerance floor.

The ordering is first-order operator splitting: `PhiEqn` and `UEqn` use the `B` field from the previous step; `B` is updated only at the end of each timestep. This is acceptable in the low-to-moderate Rm regimes targeted, where the per-step change in `B` is small.

## File layout

| File | Purpose |
| --- | --- |
| `epotMhdFoamHybrid.C` | Main time loop |
| `createFields.H` | Field declarations (`U`, `p`, `B`, `pB`, `Phi`, `phiB`, `D_B`) |
| `createControl.H` | PISO control object |
| `PhiEqn.H` | Electric-potential Poisson solve |
| `BEqn.H` | Magnetic induction equation |
| `pBEqn.H` | Helmholtz projection enforcing div(B) = 0 |
| `Make/files` | wmake target spec |
| `Make/options` | Include paths and link libraries |

## Compile

In WSL with OpenFOAM 13 sourced:

```bash
cd $WM_PROJECT_USER_DIR/applications/solvers/epotMhdFoamHybrid
wclean
wmake
```

The executable lands in `$FOAM_USER_APPBIN/epotMhdFoamHybrid`.

## Required case-file additions (vs. plain epotMhdFoam)

When migrating an existing `epotMhdFoam` case to this hybrid solver:

- **`0/pB`** — magnetic-pressure scalar field. Dimensions `[1 1 -2 0 0 -1 0]` (T·m). Suggested BCs: zero-gradient at insulating walls and inner electrode, `fixedValue uniform 0` at one electrode (anchors the projection nullspace). Initial value: `uniform 0`.
- **`system/fvSolution`** — add `pB` and `pBFinal` blocks (GAMG with low tolerance) and a `B` block (PBiCGStab / DILU). The exact tunings used in the validation campaign are in the project case directories.
- **`system/fvSchemes`** — add `div(phi,B)` and `div(phiB,U)` entries under `divSchemes`. `Gauss linear` is the safe default for both.
- **`constant/physicalProperties`** — confirm `mu` (magnetic permeability) is defined. Default vacuum value: `1.2566e-6` H/m.

## Validation status of the code

Two-stage validation  on the Corbino-disk debug mesh:

### Stage 1 — Low-Rm consistency (Rm = 0.035)

Run from a converged laminar baseline for 5 ms:

- `max |B - B_0| / B_0 = 2.57e-04` — essentially numerical noise, confirming the hybrid recovers the low-Rm limit.
- Continuity errors at `1e-11` (same level as `epotMhdFoam`).
- `max |div B| / |B|` bounded around `1.3e-04`, decreasing over the run.
- No NaN, no FPE, time loop closes cleanly.

### Stage 2 — Finite-Rm activation (Rm ~ 1.0)

Same case but with `mu` artificially scaled by 30x in `physicalProperties`, pushing Rm from 0.035 to ~1. Has been run for 50 ms:

- `max |B - B_0| / B_0 = 9.97e-03` — a 39x growth versus Stage 1.
- The 39x growth is in close agreement with the 30x Rm scaling, consistent with leading-order linear induction theory (`|B'|/B_0 ~ Rm`).
- Continuity errors remain at `1e-11`.
- `max |div B| / |B|` bounded around `5e-03`, no runaway.
- Solver remains numerically stable across the Rm range.

These results confirm the hybrid solver behaves correctly across the regimes for which it is intended.

## Sources

1. P. A. Davidson, *An Introduction to Magnetohydrodynamics*, Cambridge University Press (2001).
2. U. Müller and L. Bühler, *Magnetohydrodynamics in Channels and Containers*, Springer-Verlag (2001).
3. OpenFOAM Foundation, stock `mhdFoam` solver — the induction-equation discretisation is adapted from this source.
4. C. J. Greenshields, *OpenFOAM User Guide*, OpenFOAM Foundation (2022).

## Citation

If you use this solver in published work, please cite:

> Kocaman, A. E. (2026). *Stability of Laminar Magnetohydrodynamic Corbino Flow and Foundations for Turbulent Numerical Simulation.* Senior Project, Department of Physics, Bilkent University.

## License

Inherits from OpenFOAM (GPL v3). Distributed under the same terms.
