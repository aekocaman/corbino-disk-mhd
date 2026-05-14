# MHD Turbulence in a Corbino Disk — Senior Project

Senior physics project at Bilkent University investigating
magnetohydrodynamic (MHD) flow and stability of liquid Galinstan
in a Corbino-disk geometry under a uniform axial magnetic field.
Advisor: Aydın Cem Keser.

## Contents

- `epotMhdFoam/` — custom low-Rm electric-potential MHD solver
  (extension of OpenFOAM's stock `mhdFoam`).
- `epotMhdFoamHybrid/` — hybrid solver adding the magnetic-induction
  equation for moderate-Rm regimes. Validated against `epotMhdFoam`
  at low Rm and against linear induction theory at finite Rm.
- `postproc/` — general-purpose post-processing scripts
  (torque-closure check, Hartmann-profile fit, field diagnostics,
  probe trace analysis).
- `corbino_turbulance_debug/` — debug-mesh case files and analysis
  scripts (laminar baseline, perturbation bracket, back-EMF
  linearity check).
- `corbino_turbulance_debug_hartmann/` — Hartmann-resolved debug
  variant plus the Chebyshev-collocation Q2D linear-stability solver.
- `corbino_turbulance/` — production-mesh case files staged for
  cluster deployment (SLURM submission template, runbook,
  preflight checklist, verification demo).
- `report/` — thesis LaTeX source, headline figure script, and
  campaign log.

## Quick start

In an OpenFOAM 13 environment:

    cd epotMhdFoam       && wmake && cd ..
    cd epotMhdFoamHybrid && wmake && cd ..

See each subdirectory's README for case-specific instructions.

## Citation

Kocaman, A. E. (2026). Stability of Laminar Magnetohydrodynamic
Corbino Flow and Foundations for Turbulent Numerical Simulation.
Senior Project, Department of Physics, Bilkent University.

## License

GPL v3 (inherited from OpenFOAM).
