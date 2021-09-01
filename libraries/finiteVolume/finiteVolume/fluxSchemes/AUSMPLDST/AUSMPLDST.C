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

References
-----------
    MS. Liu (1996). A sequel to AUSM: AUSM+. J Comp Physics, 
    129, 364-382.

\*-------------------------------------------------------------------------------------------*/

#include "AUSMPLDST.H"
#include "addToRunTimeSelectionTable.H"
#include "DST.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(AUSMPLDST, 0);
addToRunTimeSelectionTable(baseFluxDST, AUSMPLDST, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AUSMPLDST::AUSMPLDST
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
            dimensionedScalar("DST::neg(", dimless, -1.0)
        )
    )        
        
{} // End of constructor

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

AUSMPLDST::~AUSMPLDST()
{
    //delete objComputeFlux;
} //end of destructor

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void AUSMPLDST::computeFluxDST()
{
        surfaceScalarField rho_own(fvc::interpolate(rho_, own_, "reconstruct(rho)"));
        surfaceScalarField rho_nei(fvc::interpolate(rho_, nei_, "reconstruct(rho)"));

        surfaceVectorField rhoU_own(fvc::interpolate(rhoU_, own_, "reconstruct(U)"));
        surfaceVectorField rhoU_nei(fvc::interpolate(rhoU_, nei_, "reconstruct(U)"));

    surfaceVectorField U_own("U_own", rhoU_own/rho_own);
        surfaceVectorField U_nei("U_nei", rhoU_nei/rho_nei);
        
        surfaceScalarField Uv_own("Uv_own", U_own & mesh_.Sf()/mesh_.magSf());
        surfaceScalarField Uv_nei("Uv_nei", U_nei & mesh_.Sf()/mesh_.magSf());

    volScalarField rPsi("rPsi", 1.0/thermo_.psi());
        surfaceScalarField rPsi_own(fvc::interpolate(rPsi, own_, "reconstruct(T)"));
        surfaceScalarField rPsi_nei(fvc::interpolate(rPsi, nei_, "reconstruct(T)"));

    volScalarField e("e", thermo_.h()-thermo_.p()/rho_); // -> stable
    //  volScalarField e("e", rhoE_/(rho_)-0.5*magSqr(rhoU_/(rho_))); -> unstable
        surfaceScalarField e_own(fvc::interpolate(e, own_, "reconstruct(T)"));
        surfaceScalarField e_nei(fvc::interpolate(e, nei_, "reconstruct(T)"));

    surfaceScalarField rhoE_own("rhoE_own",rho_own*(e_own + 0.5*magSqr(U_own)));
        surfaceScalarField rhoE_nei("rhoE_nei",rho_nei*(e_nei + 0.5*magSqr(U_nei)));

    surfaceScalarField p_own("p_own", rho_own*rPsi_own);
        surfaceScalarField p_nei("p_nei", rho_nei*rPsi_nei);

    tmp <volScalarField> H
        (
            (max(rhoE_/(rho_),dimensionedScalar("0", rhoE_.dimensions()/rho_.dimensions(), SMALL)) +
            max(thermo_.p()/(rho_),dimensionedScalar("0", thermo_.p().dimensions()/rho_.dimensions(), SMALL)))
        );
    tmp <volScalarField> a=sqrt(2.0*(thermo_.Cp()/thermo_.Cv()-1.0)/(thermo_.Cp()/thermo_.Cv()+1.0)*H);

    surfaceScalarField a_own(fvc::interpolate(a(), own_,  "reconstruct(T)"));
    surfaceScalarField a_nei(fvc::interpolate(a(), nei_,  "reconstruct(T)"));

    tmp <surfaceScalarField> a12=0.5*(a_own+a_nei);
    
        surfaceScalarField Ma_own ("Ma_own", Uv_own/a12());
        surfaceScalarField Ma_nei("Ma_nei", Uv_nei/a12());

    surfaceScalarField Ma4_own
    (
        "Ma4_own",
        0.5*(sign(1.0-mag(Ma_own))+mag(sign(1.0-mag(Ma_own))))*(0.25*sqr(Ma_own+1.0))*(1.0 - 16.0*beta_*(-0.25)*sqr(Ma_own-1.0))
        +0.5*(sign(-1.0+mag(Ma_own))+mag(sign(-1.0+mag(Ma_own))))*max(Ma_own,0.0)
    );

    surfaceScalarField P5_own
    (
        "P5_own",
        0.5*(sign(1.0-mag(Ma_own))+mag(sign(1.0-mag(Ma_own))))*0.25*sqr(Ma_own+1.0)*(2.0 - Ma_own - 16.0*alpha_*Ma_own*(-0.25)*sqr(Ma_own-1.0))
        +0.5*(sign(-1.0+mag(Ma_own))+mag(sign(-1.0+mag(Ma_own))))*DST::pos(Ma_own)
    );
    
    surfaceScalarField Ma4_nei
    (
        "Ma4_nei", 
        0.5*(sign(1.0-mag(Ma_nei))+mag(sign(1.0-mag(Ma_nei))))*(-0.25)*sqr(Ma_nei-1.0)*(1.0 + 16.0*beta_*0.25*sqr(Ma_nei+1.0))
        +0.5*(sign(-1.0+mag(Ma_nei))+mag(sign(-1.0+mag(Ma_nei))))*min(Ma_nei,0.0)
    );        

    surfaceScalarField P5_nei
    (
        "P5_nei",
        0.5*(sign(1.0-mag(Ma_nei))+mag(sign(1.0-mag(Ma_nei))))*(-0.25)*sqr(Ma_nei-1.0)*(-2.0 - Ma_nei + 16.0*alpha_*Ma_nei*0.25*sqr(Ma_nei+1.0))
        +0.5*(sign(-1.0+mag(Ma_nei))+mag(sign(-1.0+mag(Ma_nei))))*DST::neg(Ma_nei)
    );

    surfaceScalarField Ma12("Ma12", Ma4_own+Ma4_nei);

    tmp <surfaceScalarField> P12=P5_own*p_own+P5_nei*p_nei;

    Up_=DST::pos(Ma12)*U_own+DST::neg0(Ma12)*U_nei;

    phi_=Ma12*a12*(DST::pos(Ma12)*rho_own+DST::neg0(Ma12)*rho_nei)*mesh_.magSf();

    phiUp_=phi_*Up_+P12*mesh_.Sf();

    tmp <surfaceScalarField> H_own=(rhoE_own + p_own)/rho_own;
    tmp <surfaceScalarField> H_nei=(rhoE_nei + p_nei)/rho_nei;
    
    phiEp_=phi_*(DST::pos(Ma12)*H_own+DST::neg0(Ma12)*H_nei);

    calcFields();
        JacobianPtr_().computeJacobianDST();    
        
} // End of computeFluxDST() function

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
