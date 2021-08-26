/*-----------------------------------------------------------------------------------------*\
| =========                 |                                                               |
| \\      /  F ield         | rhoDST: Block-coupled solver for compressible flows           |
|  \\    /   O peration     |                                                               |
|   \\  /    A nd           | Copyright (C) 2020 Design and Simulation Tech. Inc. (DSTECH)  |
|    \\/     M anipulation  | Website:  http://www.dstechno.net/                            |
\*-----------------------------------------------------------------------------------------*/
/*
License
    rhoDST is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    rhoDST is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with rhoDST.  If not, see <http://www.gnu.org/licenses/>.

rhoDST: 
       Compressible flow solver for high speed viscous flows over 
       aeronautical vehicles. 


Author:
    Design and Simulation Technologies Inc. (DSTECH)
    http://dstechno.net/
         
      _____   _____ _______ ______ _____ _    _ 
     |  __ \ / ____|__   __|  ____/ ____| |  | |
     | |  | | (___    | |  | |__ | |    | |__| |
     | |  | |\___ \   | |  |  __|| |    |  __  |
     | |__| |____) |  | |  | |___| |____| |  | |
     |_____/|_____/   |_|  |______\_____|_|  |_|

\*---------------------------------------------------------------------------------------*/

#include "baseFluxDST.H"

namespace Foam
{

defineTypeNameAndDebug(baseFluxDST, 0);
defineRunTimeSelectionTable(baseFluxDST, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

baseFluxDST::baseFluxDST
(
    const word& type,
    const IOdictionary& dict,
    const fvMesh& mesh,
    basicThermo& thermo,
    const volScalarField& rho,
    volVectorField& U,
    const volVectorField& rhoU,
    const volScalarField& rhoE,
    compressible::turbulenceModel& turbulence,
    surfaceScalarField& phi
)
:
    JacobianPtr_
    (
        JacobianDST<tensor5>::New
        (
            dict,
            mesh,
            thermo,
            turbulence,
            rho,
            rhoU,
            rhoE
        )
    ),
    inviscid_(false),
    mesh_(mesh),
    thermo_(thermo),
    rho_(rho),
    U_(U),
    rhoU_(rhoU),
    rhoE_(rhoE),
    turbulence_(turbulence),
    phi_(phi),
    phiUp_
    (
        IOobject
        (
            "phiUp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimArea*rhoU_.dimensions()*U_.dimensions()
    ),
    phiEp_
    (
        IOobject
        (
            "phiEp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimArea*rhoE_.dimensions()*U_.dimensions()
    ),
    Up_
    (
        IOobject
        (
            "Up",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        U_.dimensions()
    ),
    R_
    (
        IOobject
        (
            "R_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector5
        (
            "R_",
            dimless,
            vector5::zero
        )
    ),
    own_
    (
        surfaceScalarField
        (
            IOobject
            (
                "own",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("own", dimless, 1.0)
        )
    ),
    nei_
    (
        surfaceScalarField
        (
            IOobject
            (
                "nei",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("neg", dimless, -1.0)
        )
    )
{} // End of constructor

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

baseFluxDST::~baseFluxDST()
{} // End of destructor

// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

surfaceScalarField& baseFluxDST::phi()
{
    return phi_;
}
surfaceVectorField& baseFluxDST::phiUp()
{
    return phiUp_;
}
surfaceScalarField& baseFluxDST::phiEp()
{
    return phiEp_;
}
surfaceVectorField& baseFluxDST::Up()
{
    return Up_;
}

JacobianDST<tensor5>& baseFluxDST::Jacobian()
{
    return JacobianPtr_();
}

vector5Field& baseFluxDST::R()
{
    return R_.internalField();
}

vector5Field& baseFluxDST::RoldTime()
{
    return R_.oldTime().internalField();
}

bool& baseFluxDST::inviscid()
{
    return inviscid_;
}

const fvMesh& baseFluxDST::mesh()
{
    return mesh_;
}

void baseFluxDST::calcFields()
{
    tmp<volScalarField> phi0(fvc::div(phi_));
    tmp<volVectorField> phi1(fvc::div(phiUp_));
    tmp<volScalarField> phi4(fvc::div(phiEp_));
    
    if(!inviscid_)
    {
        volScalarField muEff("muEff", turbulence_.muEff());
        volTensorField tauMC
        (
            "tauMC", 
            muEff*dev2(Foam::T(fvc::grad(U_)))
        );


        surfaceScalarField sigmaDotU
        (
            "sigmaDotU",
            (
                fvc::interpolate(muEff)*mesh_.magSf()*fvc::snGrad(U_)
                + (mesh_.Sf() & fvc::interpolate(tauMC))
            )
            & Up_
        );
        phi1() -= fvc::laplacian(muEff, U_) + fvc::div(tauMC);
        phi4() -= fvc::div(sigmaDotU);
    } 
    
    DST::insert(R_, phi0, 0);
    DST::insert(R_, phi1, 1) ;
    DST::insert(R_, phi4, 4);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam






