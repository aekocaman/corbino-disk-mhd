# epotMhdFoam — low-Rm electric-potential MHD solver for OpenFOAM 13

Extension of OpenFOAM Foundation's induction-only `mhdFoam` (which sits in
`legacy/electromagnetics`) to handle electrically driven liquid-metal flows
in the low magnetic-Reynolds-number limit, e.g. Corbino-disk MHD.

## Why this exists

Foundation 13's `mhdFoam` solves only the induction equation for B and uses
the resulting Maxwell stress as the body force on the fluid. That requires
B to have non-zero curl in the fluid (so J = curl(B)/mu_0 is non-zero), which
is true for some Hartmann-channel setups but **not** for a Corbino disk with
a uniform applied B and current driven by an electrode potential V. In our
configuration the original solver produced no body force at all (confirmed empirically: an unfortunate 8 wall-clock hours yielded
U ~ 1e-13 m/s).

This solver solves a Poisson equation for the electric potential Phi,
constructs the current density J = sigma * (-grad(Phi) + U x B), and uses
J x B directly as the Lorentz force in the momentum equation. B is treated
as fixed at the imposed value (low-Rm: induction is negligible).

Validity range: Rm = mu_0 * sigma * U * L << 1. 

## Files

```
epotMhdFoam/
    Make/
        files               compile target name
        options             include paths and link libraries
    epotMhdFoam.C           main solver loop
    createControl.H         pisoControl instantiation
    createFields.H          field & property registration (adds Phi)
    PhiEqn.H                electric-potential Poisson solve
```

## Installation in WSL

```bash
# 1. Place the source in the user solver directory
mkdir -p 
cp -r 

# 2. Source the OpenFOAM 13 environment if not already in your bashrc
source /opt/openfoam13/etc/bashrc

# 3. Compile
cd ~/OpenFOAM/reybo78-13/applications/solvers/epotMhdFoam
wmake


```


## Run the case

The existing `case_360_laminar/` 's `0/Phi` and `system/fvSolution`  are essentially compatible.
 Two changes are needed:

1. **`system/fvSchemes`**: add divergence and laplacian schemes for the new
   terms if not already present:

       laplacianSchemes
       {
           ...
           laplacian(sigma,Phi)  Gauss linear corrected;
       }

       divSchemes
       {
           ...
           div(sigma*(U^B))      Gauss linear;
       }

2. **Run the new solver instead of `mhdFoam`**:

       cd ~/OpenFOAM/reybo78-13/run/case_360_laminar
       #  to remove old snapshots 
       foamListTimes -rm
       # restart from t=0
       nohup epotMhdFoam > log.epotMhdFoam 2>&1 &
       echo $! > epotMhdFoam.pid


