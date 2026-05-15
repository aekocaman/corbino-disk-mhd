// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
// Electric-potential Poisson equation under low-Rm approximation.
//
// Charge conservation div(J) = 0 with Ohm's law J = sigma*(-grad(Phi) + U x B)
// gives, for spatially uniform sigma:
//
//      div( sigma * grad(Phi) ) = div( sigma * (U x B) )
//
// Boundary conditions are set in 0/Phi
//
// Solver settings (tolerances, preconditioner) come from system/fvSolution
// under the "Phi" entry.
//
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

{
    fvScalarMatrix PhiEqn
    (
        fvm::laplacian(sigma, Phi)
     ==
        fvc::div(sigma*(U ^ B))
    );

    PhiEqn.solve();
}
