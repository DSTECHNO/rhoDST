Info<< "Reading thermophysical properties\n" << endl;

autoPtr<basicPsiThermo> thermo
(
    basicPsiThermo::New(mesh)
);

volScalarField& h = thermo->h();
volScalarField& p = thermo->p();

const volScalarField& T = thermo->T();
const volScalarField& psi = thermo->psi();
const volScalarField& mu = thermo->mu();

bool inviscid(true);
if (max(mu.internalField()) > 0.0)
{
    inviscid = false;
}

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
    thermo->rho()
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

volVectorField rhoU
(
    IOobject
    (
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U
);

volScalarField rhoE
(
    IOobject
    (
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*(h*thermo->Cv()/thermo->Cp() + 0.5*magSqr(U))
);

#include "compressibleCreatePhi.H"

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);    

Info<< "Reading thermophysicalProperties\n" << endl;

// Pr defined as a separate constant to enable calculation of k, currently
// inaccessible through thermo
IOdictionary thermophysicalProperties
(
    IOobject
    (
        "thermophysicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar Pr = dimensionedScalar("Pr", dimless, 1.0);

if (thermophysicalProperties.found("Pr"))
{
    Pr = dimensionedScalar(thermophysicalProperties.lookup("Pr"));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<Foam::baseFluxDST> rhoDSTFluxPtr;

rhoDSTFluxPtr = 
(
    Foam::baseFluxDST::New
    (
        rhoDSTDict,
        mesh,
        thermo(),
        rho,
        U,
        rhoU,
        rhoE,
        turbulence(),
        phi
    )
);

Foam::baseFluxDST& rhoDSTFlux = rhoDSTFluxPtr(); 

runTime.initialize
(
    thermo(),
    rhoDSTDict,
    turbulence(),
    U
);

autoPtr<iteration> rhoDSTPtr
(
    iteration::New
    (
        rhoDSTDict,
        runTime,
        thermo(),
        U,
        turbulence(),
        rhoDSTFlux
    )
);
iteration& rhoDST = rhoDSTPtr();

volVector5FieldDST& W = rhoDST.getW();
W.set(rho, rhoU, rhoE);
