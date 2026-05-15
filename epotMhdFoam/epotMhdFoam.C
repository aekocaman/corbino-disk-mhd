/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | OpenFOAM 13 user-extension
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    epotMhdFoam

Description
    Low-magnetic-Reynolds-number incompressible MHD solver for an
    electrically conducting fluid driven by an applied electrostatic
    potential at the boundaries (e.g. the Corbino-disk geometry).

    Suitable when Rm = mu_0 * sigma * U * L << 1, so that induction of B by
    fluid motion is negligible and B is treated as fixed.  The current density is

        J = sigma * (-grad(Phi) + U x B)

    and the electric potential Phi satisfies charge conservation:

        div(sigma * grad(Phi)) = div(sigma * (U x B))

    The Lorentz body force F = J x B is decomposed as

        F = sigma * ( -grad(Phi) x B + (U.B) B - |B|^2 U )

    The dissipative -sigma * |B|^2 * U term (Hartmann damping) is treated
    implicitly via fvm::Sp; this is essential for stability at high
    Hartmann number, where the explicit timestep would otherwise be
    bounded by the magnetic-damping timescale tau_H ~ rho/(sigma B^2).

    It is an modified version of OpenFOAM Foundation's induction-only mhdFoam (legacy/
    electromagnetics).  The induction equation, the magnetic-pressure
    Poisson step, and the BPISO loop are removed; a Phi Poisson equation
    is added; and the Maxwell-stress force terms in UEqn are replaced by
    the J x B Lorentz force.

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

    Info<< nl << "Starting time loop" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        #include "CourantNo.H"

        // ---------------------------------------------------------------- //
        //  Adaptive Dt control (Courant-based)
        //
        //  mhdFoam's loop has no setDeltaT.H, so we inline the standard
        //  controller here.  Reading the controls each step lets the user
        //  retune adjustTimeStep / maxCo / maxDeltaT mid-run via
        //  runTimeModifiable yes.
        // ---------------------------------------------------------------- //
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
        //  Electric-potential Poisson solve
        // ---------------------------------------------------------------- //
        #include "PhiEqn.H"

        // ---------------------------------------------------------------- //
        //  Diagnostic current density and Lorentz-force decomposition
        // ---------------------------------------------------------------- //
        // Auto-write so postProcess can integrate J for torque diagnostics.
        volVectorField J
        (
            IOobject
            (
                "J",
                runTime.name(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sigma*(-fvc::grad(Phi) + (U ^ B))
        );

        // -sigma * |B|^2 * U  (implicit Hartmann damping coefficient)
        volScalarField sigmaBsqByRho
        (
            "sigmaBsqByRho",
            sigma*magSqr(B)/rho
        );

        // sigma * (-grad(Phi) x B) / rho       — driving force from Phi
        volVectorField F_drive
        (
            "F_drive",
            (-sigma*(fvc::grad(Phi) ^ B))/rho
        );

        // sigma * (U.B) * B / rho              — small for axial B
        volVectorField F_UBdotB
        (
            "F_UBdotB",
            (sigma*(U & B)*B)/rho
        );

        // ---------------------------------------------------------------- //
        //  Momentum predictor with implicit Hartmann damping
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

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
