/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | OpenFOAM 13 user-extension
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    epotMhdFoamHybrid

Description
    Hybrid epotMhdFoam-mhdFoam solver: the low-Rm electric-potential
    formulation of epotMhdFoam combined with the induction equation of
    stock mhdFoam.  Suitable for MHD flows ranging from low to moderate
    magnetic Reynolds number (Rm <~ O(1)).

    Equations solved (per timestep):

      1. Electric-potential Poisson (uses current total B):
             div(sigma * grad(Phi)) = div(sigma * (U x B))

      2. Lorentz body force decomposition (uses current total B):
             F_L / rho = sigma/rho * [-grad(Phi) x B + (U.B) B - |B|^2 U]

      3. Momentum predictor + PISO loop for U and p (Hartmann damping
         -sigma|B|^2 U / rho treated implicitly via fvm::Sp).

      4. Induction equation for B:
             dt(B) + div(phi,B) - div(phiB,U) - laplacian(DB,B) = 0
         with DB = 1/(mu * sigma).

      5. Solenoidal projection of B via magnetic-pressure pB:
             laplacian(pB) = div(phiB)
             phiB := phiB - flux(grad(pB)),    B := B - grad(pB)

    The ordering is operator-splitting at first order: PhiEqn and UEqn
    use the B from the previous step; B is updated only at the end.
    Acceptable at low Rm where B' contributions per step are small.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "pisoControl.H"
#include "pressureReference.H"
#include "constrainPressure.H"
#include "constrainHbyA.H"
#include "fvcDdt.H"
#include "fvcGrad.H"
#include "fvcFlux.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "Starting time loop (hybrid solver: low-Rm + induction)" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        #include "CourantNo.H"

        // Adaptive Dt (Courant-based), copied from epotMhdFoam
        {
            const bool adjustTimeStep =
                runTime.controlDict().lookupOrDefault("adjustTimeStep", false);
            const scalar maxCo =
                runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);
            const scalar maxDeltaT =
                runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

            if (adjustTimeStep)
            {
                const scalar maxDeltaTFact = maxCo/(CoNum + SMALL);
                const scalar deltaTFact =
                    min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);
                runTime.setDeltaT
                (
                    min(deltaTFact*runTime.deltaTValue(), maxDeltaT)
                );
                Info<< "deltaT = " << runTime.deltaTValue() << endl;
            }
        }

        // ---------------------------------------------------------------- //
        //  1.  Electric-potential Poisson  (uses current B)
        // ---------------------------------------------------------------- //
        #include "PhiEqn.H"

        // ---------------------------------------------------------------- //
        //  2.  Diagnostic J + Lorentz-force decomposition (uses current B)
        // ---------------------------------------------------------------- //
        volVectorField J
        (
            IOobject("J", runTime.name(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
            sigma*(-fvc::grad(Phi) + (U ^ B))
        );

        volScalarField sigmaBsqByRho("sigmaBsqByRho", sigma*magSqr(B)/rho);
        volVectorField F_drive("F_drive", (-sigma*(fvc::grad(Phi) ^ B))/rho);
        volVectorField F_UBdotB("F_UBdotB", (sigma*(U & B)*B)/rho);

        // ---------------------------------------------------------------- //
        //  3.  Momentum predictor with implicit Hartmann damping
        // ---------------------------------------------------------------- //
        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
          + fvm::Sp(sigmaBsqByRho, U)
         ==
            F_drive + F_UBdotB
        );

        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
        }

        // ---------------------------------------------------------------- //
        //  PISO loop (pressure correction and U update)
        // ---------------------------------------------------------------- //
        while (piso.correct())
        {
            volScalarField rAU(1.0/UEqn.A());
            surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::flux(HbyA)
              + rAUf*fvc::ddtCorr(U, phi)
            );

            constrainPressure(p, U, phiHbyA, rAUf);

            while (piso.correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAUf, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference
                (
                    pressureReference.refCell(),
                    pressureReference.refValue()
                );

                pEqn.solve();

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
        }

        // ---------------------------------------------------------------- //
        //  4.  Induction equation for B (NEW)
        // ---------------------------------------------------------------- //
        #include "BEqn.H"

        // ---------------------------------------------------------------- //
        //  5.  Solenoidal projection of B via magnetic pressure pB (NEW)
        // ---------------------------------------------------------------- //
        #include "pBEqn.H"

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
