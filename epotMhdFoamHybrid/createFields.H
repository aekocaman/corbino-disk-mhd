// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
// Field initialisation for epotMhdFoamHybrid.
//
// Same as epotMhdFoam createFields.H, plus:
//   - pB    : magnetic-pressure scalar (for solenoidal projection of B)
//   - phiB  : face flux of B (computed from B; used in BEqn and pBEqn)
//   - DB    : magnetic diffusivity 1/(mu * sigma)
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

dimensionedScalar rho("rho", dimDensity, physicalProperties);
dimensionedScalar nu("nu", dimKinematicViscosity, physicalProperties);
dimensionedScalar mu("mu", dimensionSet(1, 1, -2, 0, 0, -2, 0), physicalProperties);
dimensionedScalar sigma("sigma", dimensionSet(-1, -3, 3, 0, 0, 2, 0), physicalProperties);

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject("p", runTime.name(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject("U", runTime.name(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);

#include "createPhi.H"

Info<< "Reading field B (now evolved by induction equation)\n" << endl;
volVectorField B
(
    IOobject("B", runTime.name(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);

Info<< "Reading field pB (magnetic-pressure scalar)\n" << endl;
volScalarField pB
(
    IOobject("pB", runTime.name(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);

Info<< "Reading field Phi (electric potential)\n" << endl;
volScalarField Phi
(
    IOobject("Phi", runTime.name(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);

// Magnetic face flux of B (recomputed each step from B)
surfaceScalarField phiB
(
    IOobject("phiB", runTime.name(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    fvc::flux(B)
);

// Magnetic diffusivity:  DB = 1 / (mu * sigma),  units [m^2/s]
dimensionedScalar DB = 1.0/(mu*sigma);
DB.name() = "DB";

pressureReference pressureReference(p, piso.dict());

mesh.schemes().setFluxRequired(p.name());
mesh.schemes().setFluxRequired(Phi.name());
mesh.schemes().setFluxRequired(pB.name());
