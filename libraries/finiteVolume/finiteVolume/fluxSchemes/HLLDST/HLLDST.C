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

\*---------------------------------------------------------------------------------------*/

#include "HLLDST.H"
#include "addToRunTimeSelectionTable.H"
#include "DST.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(HLLDST, 0);
addToRunTimeSelectionTable(baseFluxDST, HLLDST, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

HLLDST::HLLDST
(
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
    )        
{} // End of constructor

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

HLLDST::~HLLDST()
{} //End of destructor

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void HLLDST::computeFluxDST()
{
    surfaceScalarField rho_own(fvc::interpolate(rho_, own_, "reconstruct(rho)"));
        surfaceScalarField rho_nei(fvc::interpolate(rho_, nei_, "reconstruct(rho)"));    

    surfaceVectorField rhoU_own(fvc::interpolate(rhoU_, own_, "reconstruct(U)"));
        surfaceVectorField rhoU_nei(fvc::interpolate(rhoU_, nei_, "reconstruct(U)"));

    surfaceVectorField U_own("U_own", rhoU_own/rho_own);
        surfaceVectorField U_nei("U_nei", rhoU_nei/rho_nei);

        volScalarField e("e", thermo_.h()-thermo_.p()/rho_); // -> stable
    //  volScalarField e("e", rhoE_/(rho_)-0.5*magSqr(rhoU_/(rho_))); -> unstable
        surfaceScalarField e_own(fvc::interpolate(e, own_, "reconstruct(T)"));
        surfaceScalarField e_nei(fvc::interpolate(e, nei_, "reconstruct(T)"));

    surfaceScalarField E_own("E_own",e_own + 0.5*magSqr(U_own));
        surfaceScalarField E_nei("E_nei",e_nei + 0.5*magSqr(U_nei));

    volScalarField rPsi("rPsi", 1.0/thermo_.psi());
        surfaceScalarField rPsi_own(fvc::interpolate(rPsi, own_, "reconstruct(T)"));
        surfaceScalarField rPsi_nei(fvc::interpolate(rPsi, nei_, "reconstruct(T)"));

    surfaceScalarField p_own("p_own", rho_own*rPsi_own);
        surfaceScalarField p_nei("p_nei", rho_nei*rPsi_nei);

    surfaceScalarField H_own(E_own+p_own/rho_own); 
    surfaceScalarField H_nei(E_nei+p_nei/rho_nei);         

    //volScalarField a("a", sqrt(thermo->Cp()/thermo->Cv()*rPsi));

    tmp <volScalarField> H
        (
            (max(rhoE_/(rho_),dimensionedScalar("0", rhoE_.dimensions()/rho_.dimensions(), SMALL)) +
            max(thermo_.p()/(rho_),dimensionedScalar("0", thermo_.p().dimensions()/rho_.dimensions(), SMALL)))
        );

    volScalarField a=sqrt(2.0*(thermo_.Cp()/thermo_.Cv()-1.0)/(thermo_.Cp()/thermo_.Cv()+1.0)*H);

        surfaceScalarField a_own
        (
            "a_own",
            fvc::interpolate(a, own_, "reconstruct(T)")
        );
        surfaceScalarField a_nei
        (
            "a_nei",
            fvc::interpolate(a, nei_, "reconstruct(T)")
        );

    //Roe averaging    
    dimensionedScalar rho0("rho0", dimDensity, VSMALL);
    tmp <surfaceScalarField> coefR=sqrt(max(rho0,rho_nei)/max(rho0,rho_own));

    tmp <surfaceScalarField> HTilde
        (
                (coefR()*H_nei+H_own)/(coefR()+1.0)
        );

    tmp <surfaceVectorField> UTilde
    ( 
        (coefR()*U_nei+U_own)/(coefR()+1.0)
    );

    surfaceScalarField aTilde
    (
        sqrt((fvc::interpolate(thermo_.Cp()/thermo_.Cv())-1.0)*(HTilde()-0.5*magSqr(UTilde())))
    );

    surfaceScalarField Uv_own=U_own & mesh_.Sf()/mesh_.magSf();
    surfaceScalarField Uv_nei=U_nei & mesh_.Sf()/mesh_.magSf();
    surfaceScalarField UvTilde=UTilde() & mesh_.Sf()/mesh_.magSf();

    surfaceScalarField S_own("S_own", min(Uv_own-a_own, UvTilde-aTilde));
    surfaceScalarField S_nei("S_nei", max(Uv_nei+a_nei, UvTilde+aTilde));

    // pos and neg functions are defined in a different way in Extend and also pos0 and neg0 functions are not available in Extend. 
    // Thus, DST primitive functions are used here
    surfaceScalarField cf_own("cf_own",DST::pos0(S_own));             // S_own>=0 ? 1 : 0
    surfaceScalarField cf_m("cf_m", DST::neg(S_own)*DST::pos0(S_nei));     // S_own<0 and S_nei>=0 ? 1 : 0                              
        surfaceScalarField cf_nei("cf_nei",DST::neg(S_nei)*DST::neg(S_own));    // S_own<0 and S_nei<0 ? 1 : 0  

    surfaceVectorField rhoUPhi_own ("rhoUPhi_own", rho_own*U_own*Uv_own+p_own*mesh_.Sf()/mesh_.magSf());
    surfaceVectorField rhoUPhi_nei ("rhoUPhi_nei", rho_nei*U_nei*Uv_nei+p_nei*mesh_.Sf()/mesh_.magSf());

    phi_ = 
        (
        cf_own*rho_own*Uv_own +
              cf_m*(S_nei*rho_own*Uv_own-S_own*rho_nei*Uv_nei+S_own*S_nei*(rho_nei-rho_own))/(S_nei-S_own) + 
              cf_nei*rho_nei*Uv_nei
    )*mesh_.magSf();

        phiUp_ = 
    (
        cf_own*rhoUPhi_own +
        cf_m*(S_nei*rhoUPhi_own-S_own*rhoUPhi_nei+S_own*S_nei*(rho_nei*U_nei-rho_own*U_own))/(S_nei-S_own) + 
        cf_nei*rhoUPhi_nei
    )*mesh_.magSf(); 

    phiEp_ =
    (
        cf_own*rho_own*Uv_own*H_own +
        cf_m*(S_nei*rho_own*Uv_own*H_own-S_own*rho_nei*Uv_nei*H_nei+S_own*S_nei*(rho_nei*E_nei-rho_own*E_own))/(S_nei-S_own) +
        cf_nei*rho_nei*Uv_nei*H_nei
    )*mesh_.magSf();

    Up_ = cf_own*U_own +
         cf_m*(S_nei*U_own-S_own*U_nei+Uv_own*U_own-Uv_nei*U_nei)/(S_nei-S_own) + 
         cf_nei*U_nei;

    calcFields();
        JacobianPtr_().computeJacobianDST();  

} // End of computeFluxDST() function  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam




    

 
