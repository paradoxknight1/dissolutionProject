Info<< "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar molecularDiffusivity
(
    "D",
    dimensionSet(0, 2, -1, 0, 0, 0, 0),
    transportProperties
);

dimensionedScalar nu
(
    "nu",
    dimensionSet(0, 2, -1, 0, 0, 0, 0),
    transportProperties
);

dimensionedScalar rhoCarrier
(
    "rhoCarrier",
    dimensionSet(1, -3, 0, 0, 0, 0, 0),
    transportProperties
);

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
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
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field C\n" << endl;
volScalarField C
(
    IOobject
    (
        "C",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Creating dissolution source field S_C\n" << endl;
volScalarField S_C
(
    IOobject
    (
        "S_C",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zeroSource",
        dimensionSet(1, -3, -1, 0, 0, 0, 0),
        0
    )
);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    rhoCarrier
);

volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    rhoCarrier*nu
);

#include "readGravitationalAcceleration.H"

#include "createPhi.H"

pressureReference pressureReference(p, piso.dict());

mesh.setFluxRequired(p.name());
mesh.setFluxRequired(C.name());

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"

Info<< "Creating user-space dissolution cloud\n" << endl;
DissolutionCloud dissolutionCloud
(
    "dissolutionCloud",
    rho,
    U,
    mu,
    g,
    S_C
);
