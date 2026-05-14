// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
// Induction equation:
//      dt(B) = curl(U x B) + DB * laplacian(B)
//
// Discretised form (incompressible, solenoidal U and B):
//      dt(B) + div(phi, B) - div(phiB, U) - laplacian(DB, B) = 0
//
// where:
//   phi     = volumetric face flux of U  (from createPhi.H, updated by PISO)
//   phiB    = face flux of B            (recomputed below from current B)
//   DB      = 1/(mu * sigma)            (magnetic diffusivity)
//
// The (B.grad)U field-line stretching term is provided by -div(phiB, U).
// Diffusion is implicit; advection of B by U is implicit; "advection of U
// by B" (the field-line stretching) is explicit.
//
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Refresh phiB from the current B field before solving
phiB = fvc::flux(B);

fvVectorMatrix BEqn
(
    fvm::ddt(B)
  + fvm::div(phi, B)
  - fvc::div(phiB, U)
  - fvm::laplacian(DB, B)
);

BEqn.solve();

// phiB now needs re-projection in pBEqn.H to remove discretisation
// errors that have made it non-solenoidal during the BEqn solve.
phiB = fvc::flux(B);
