#   The abandoned version of the System

This case predates the custom epotMhdFoam solver. It was set up
using OpenFOAM's stock `mhdFoam`, which evolves the magnetic field
via the induction equation but only generates a Lorentz force from
the magnetic-field gradient (Maxwell stress). For a uniform axial
B = 0.5 T this gradient is zero, so the fluid remains essentially
at rest (|U| ~ 1e-13 m/s after an 8-hour run).

This observation motivated the development of `epotMhdFoam` (you can check), which uses the electric-potential formulation
and produces the correct Lorentz body force even for uniform B.

Do not use it as a starting point for new work. I leave it for the curious ones. See
`corbino_turbulance_debug/` for the active debug case.
