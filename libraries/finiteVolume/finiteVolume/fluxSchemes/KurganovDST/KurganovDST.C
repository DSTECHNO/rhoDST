/*-----------------------------------------------------------------------------------------*\
| =========                 |                                                               |
| \\      /  F ield         | rhoDST: Block-coupled solver for compressible flows           |
|  \\    /   O peration     |                                                               |
|   \\  /    A nd           | Copyright (C) 2020 Design and Simulation Tech. Inc. (DSTECH)  |
|    \\/     M anipulation  | Website:  http://www.dstechno.net/                            |
\*-----------------------------------------------------------------------------------------*/
/*
License
    This file is derivative work of rhoDST.

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

\*--------------------------------------------------------------------------------------------*/

#include "KurganovDST.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(KurganovDST, 0);
addToRunTimeSelectionTable(baseFluxDST, KurganovDST, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

KurganovDST::KurganovDST
(
    const IOdictDST& dict,
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
    baseFluxDST
    (
        typeName, 
        dict, 
        mesh,
        thermo, 
        rho, 
        U, 
        rhoU, 
        rhoE, 
        turbulence,
        phi
    ),
    v_zero("v_zero", dimVolume/dimTime, 0.0)
{} // End of constructor

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

KurganovDST::~KurganovDST()
{
} // End of destructor

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void KurganovDST::computeFluxDST()
{
    surfaceScalarField rho_own(fvc::interpolate(rho_, own_, "reconstruct(rho)"));
    surfaceScalarField rho_nei(fvc::interpolate(rho_, nei_, "reconstruct(rho)"));

    surfaceVectorField rhoU_own(fvc::interpolate(rhoU_, own_, "reconstruct(U)"));
    surfaceVectorField rhoU_nei(fvc::interpolate(rhoU_, nei_, "reconstruct(U)"));

    volScalarField rPsi("rPsi", 1.0/thermo_.psi());
    surfaceScalarField rPsi_own(fvc::interpolate(rPsi, own_, "reconstruct(T)"));
    surfaceScalarField rPsi_nei(fvc::interpolate(rPsi, nei_, "reconstruct(T)"));

    volScalarField e("e", thermo_.h()-thermo_.p()/rho_); // -> stable
//  volScalarField e("e", rhoE_/(rho_)-0.5*magSqr(rhoU_/(rho_)));// -> unstable
    surfaceScalarField e_own(fvc::interpolate(e, own_, "reconstruct(T)"));
    surfaceScalarField e_nei(fvc::interpolate(e, nei_, "reconstruct(T)"));

    surfaceVectorField U_own("U_own", rhoU_own/rho_own);
    surfaceVectorField U_nei("U_nei", rhoU_nei/rho_nei);

    surfaceScalarField rhoE_own("rhoE_own",rho_own*(e_own + 0.5*magSqr(U_own)));
    surfaceScalarField rhoE_nei("rhoE_nei",rho_nei*(e_nei + 0.5*magSqr(U_nei)));

    surfaceScalarField p_own("p_own", rho_own*rPsi_own);
    surfaceScalarField p_nei("p_nei", rho_nei*rPsi_nei);

    surfaceScalarField phiv_own("phiv_own", U_own & mesh_.Sf());
    surfaceScalarField phiv_nei("phiv_nei", U_nei & mesh_.Sf());

    volScalarField c("c", sqrt(thermo_.Cp()/thermo_.Cv()*rPsi));
    surfaceScalarField cSf_own
    (
        "cSf_own",
        fvc::interpolate(c, own_, "reconstruct(T)")*mesh_.magSf()
    );
    surfaceScalarField cSf_nei
    (
        "cSf_nei",
        fvc::interpolate(c, nei_, "reconstruct(T)")*mesh_.magSf()
    );

    surfaceScalarField ap
    (
        "ap",
        max(max(phiv_own + cSf_own, phiv_nei + cSf_nei), v_zero)
    );
    surfaceScalarField am
    (
        "am",
        min(min(phiv_own - cSf_own, phiv_nei - cSf_nei), v_zero)
    );

    surfaceScalarField a_own("a_own", ap/(ap - am));

    surfaceScalarField aSf("aSf", am*a_own);

    surfaceScalarField a_nei("a_nei", 1.0 - a_own);

    phiv_own *= a_own;
    phiv_nei *= a_nei;

    surfaceScalarField aphiv_own("aphiv_own", phiv_own - aSf);
    surfaceScalarField aphiv_nei("aphiv_nei", phiv_nei + aSf);

    phi_ = aphiv_own*rho_own + aphiv_nei*rho_nei;

    phiUp_=(aphiv_own*rhoU_own + aphiv_nei*rhoU_nei)+ (a_own*p_own + a_nei*p_nei)*mesh_.Sf();

    phiEp_=aphiv_own*(rhoE_own + p_own)+ aphiv_nei*(rhoE_nei + p_nei)+ aSf*p_own - aSf*p_nei;

    Up_=a_own*U_own + a_nei*U_nei;
    
    calcFields();
    JacobianPtr_().computeJacobianDST();
} // End of computeFluxDST() function

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
