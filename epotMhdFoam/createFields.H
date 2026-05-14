// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
// Field initialisation for epotMhdFoam.
//
// Reads all transport properties (nu, rho, sigma, mu) from
// constant/transportProperties, then registers p, U, B, Phi as volume fields
// and creates the volumetric face flux phi from U.
//
// Differences from mhdFoam's createFields.H:
//   - adds the electric potential Phi (volScalarField)
//   - drops DBU (1/(mu*rho)) and the related magnetic-pressure pB and phiB
//     since the induction equation is not solved in the low-Rm limit
//   - registers Phi as flux-required so the pressure projection on Phi is
//     consistent with OpenFOAM's PISO machinery
//
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Info<< "Reading physicalProperties\n" << endl;

IOdictionary physicalProperties
(
    IOobject
    (
        "physicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar rho
(
    "rho",
    dimDensity,
    physicalProperties
);

dimensionedScalar nu
(
    "nu",
    dimKinematicViscosity,
    physicalProperties
);

// mu: magnetic permeability, kept for completeness (unused in low-Rm UEqn)
dimensionedScalar mu
(
    "mu",
    dimensionSet(1, 1, -2, 0, 0, -2, 0),
    physicalProperties
);

// sigma: electrical conductivity, [S/m] = [A^2 s^3 / (kg m^3)]
dimensionedScalar sigma
(
    "sigma",
    dimensionSet(-1, -3, 3, 0, 0, 2, 0),
    physicalProperties
);

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

Info<< "Reading field B (applied magnetic field, fixed in low-Rm limit)\n"
    << endl;
volVectorField B
(
    IOobject
    (
        "B",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Phi (electric potential)\n" << endl;
volScalarField Phi
(
    IOobject
    (
        "Phi",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

pressureReference pressureReference(p, piso.dict());

mesh.schemes().setFluxRequired(p.name());
mesh.schemes().setFluxRequired(Phi.name());
