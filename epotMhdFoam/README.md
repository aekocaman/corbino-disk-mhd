# epotMhdFoam — low-Rm electric-potential MHD solver for OpenFOAM 13

Extension of OpenFOAM Foundation's induction-only `mhdFoam` (which sits in
`legacy/electromagnetics`) to handle electrically driven liquid-metal flows
in the low magnetic-Reynolds-number limit, e.g. Corbino-disk MHD.

## Why this exists

Foundation 13's `mhdFoam` solves only the induction equation for B and uses
the resulting Maxwell stress as the body force on the fluid. That requires
B to have non-zero curl in the fluid (so J = curl(B)/mu_0 is non-zero), which
is true for some Hartmann-channel setups but **not** for a Corbino disk with
a uniform applied B and current driven by an electrode potential V. In that
configuration the original solver produces no body force at all, no flow,
and nothing useful (confirmed empirically: 8 wall-clock hours yielded
U ~ 1e-13 m/s).

This solver solves a Poisson equation for the electric potential Phi,
constructs the current density J = sigma * (-grad(Phi) + U x B), and uses
J x B directly as the Lorentz force in the momentum equation. B is treated
as fixed at the imposed value (low-Rm: induction is negligible).

Validity range: Rm = mu_0 * sigma * U * L << 1. For Galinstan in your
geometry: Rm ~ 4*pi*1e-7 * 3.46e6 * 0.1 * 0.08 ~ 0.035, well within the
low-Rm regime.

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

## Install and build (in WSL)

```bash
# 1. Place the source in the user solver directory
mkdir -p ~/OpenFOAM/reybo78-13/applications/solvers
cp -r /mnt/c/Users/user/Documents/Claude/Projects/Turbulance\ Instabilities\ of\ a\ conducting\ fluid\ in\ 3-D\ Corbino\ Geomtry/epotMhdFoam \
      ~/OpenFOAM/reybo78-13/applications/solvers/

# 2. Source the OpenFOAM 13 environment if not already in your bashrc
source /opt/openfoam13/etc/bashrc

# 3. Compile
cd ~/OpenFOAM/reybo78-13/applications/solvers/epotMhdFoam
wmake

# 4. Confirm the binary
which epotMhdFoam
# expected: ~/OpenFOAM/reybo78-13/platforms/linux64GccDPInt32Opt/bin/epotMhdFoam
```

If `wmake` fails, the typical first culprits and fixes:

- **Header not found**: a Foundation API header path is different from the
  guess in `Make/options`. Add the missing include directory to `EXE_INC`.
- **Undefined symbol**: a runtime library is missing. Add it to `EXE_LIBS`
  with `-l<libname>`.
- **`fvm::Sp` overload**: if the compiler complains about
  `fvm::Sp(volScalarField, volVectorField)`, replace it with the equivalent
  explicit form `-= sigmaBsqByRho * U` and remove the `fvm::Sp` line.
  This sacrifices the implicit-damping benefit but compiles trivially.

## Run the case

The existing `case_360_laminar/` is essentially compatible — it already has
`0/Phi` with the right BCs, `system/fvSolution` already has a Phi solver
block. Two changes are needed:

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
       # remove old snapshots (they were physics-meaningless anyway)
       foamListTimes -rm
       # restart from t=0
       nohup epotMhdFoam > log.epotMhdFoam 2>&1 &
       echo $! > epotMhdFoam.pid

## Sanity checks once running

Within the first few timesteps:

```bash
grep "Solving for Phi" log.epotMhdFoam | tail -5
# Should now appear -- proves the Phi equation is being integrated.

grep "Solving for U" log.epotMhdFoam | tail -5
# Initial residuals should grow into ~1e-3 to 1e-2 range
# (vs the 1e-5 noise floor of the old mhdFoam runs).

cat postProcessing/probes/0/U | head -30
# U at the midplane probes should reach ~mm/s within ms.
```

Once steady state is reached, the closure check:

```
EM torque on inner electrode  +  EM torque on outer electrode
                          =  viscous torque on insulating plates
```

(both signs flipped) should hold to within ~1%.
