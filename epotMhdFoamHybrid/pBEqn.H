// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
// Solenoidal projection of B:
//      laplacian(pB) = div(phiB)            -> solve for pB
//      phiB := phiB - flux(grad(pB))        -> face flux made divergence-free
//      B    := B    - grad(pB)              -> cell B updated to match
//
// After this step, |div(B)|/|B| should drop to machine precision in the
// interior.  Boundary B values are left untouched by their fixed-value BCs.
//
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

{
    fvScalarMatrix pBEqn
    (
        fvm::laplacian(pB) == fvc::div(phiB)
    );

    pBEqn.solve();

    // Project the face flux to be divergence-free
    phiB -= pBEqn.flux();

    // Update cell-centred B consistently
    B -= fvc::grad(pB);
    B.correctBoundaryConditions();

    // Diagnostic: max(|div B|/|B|) over all cells
    const dimensionedScalar smallB("smallB", B.dimensions(), SMALL);
    const volScalarField errField
    (
        mag(fvc::div(phiB))/(mag(B) + smallB)
    );
    Info<< "Magnetic field error: max(|div B|/|B|) = "
        << max(errField).value() << endl;
}
