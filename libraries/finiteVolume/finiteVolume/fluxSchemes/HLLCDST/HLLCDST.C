/*-----------------------------------------------------------------------------------------*\
| =========                 |                                                               |
| \\      /  F ield         | rhoDST: Block-coupled solver for compressible flows           |
|  \\    /   O peration     |                                                               |
|   \\  /    A nd           | Copyright (C) 2020 Design and Simulation Tech. Inc. (DSTECH)  |
|    \\/     M anipulation  | Website:  http://www.dstechno.net/                            |
\*-----------------------------------------------------------------------------------------*/
/*
License
    This file is derivative work of OpenrhoDST.

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
    Eleuterio F. Toro, Michael Spruce, and William Speares. Restoration of the
     contact surface in the HLL-Riemann solver. Shock waves, 4(1):25–34, 1994.

\*---------------------------------------------------------------------------------------*/


#include "HLLCDST.H"
#include "addToRunTimeSelectionTable.H"
#include "DST.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(HLLCDST, 0);
addToRunTimeSelectionTable(baseFluxDST, HLLCDST, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

HLLCDST::HLLCDST
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

    own_(surfaceScalarField
        (
            IOobject
            (
                "own",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("own", dimless, 1.0)
        )),

        nei_(surfaceScalarField
        (
            IOobject
            (
                "nei",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("neg", dimless, -1.0)
        ))        
        
        
{
    directWaveSpeed=mesh_.lookupObject<IOdictionary>("rhoDSTDict").lookupOrDefault("directWaveSpeed", true);        

    if (directWaveSpeed)
    {
        Info << "-----> Direct wave speed estimates" << endl;
    }
    else
        Info << "-----> pressure-based wave speed estimates" << endl;


} // End of constructor

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

HLLCDST::~HLLCDST()
{} //End of destructor

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void HLLCDST::computeFluxDST()
{
    surfaceScalarField rho_own(fvc::interpolate(rho_, own_, "reconstruct(rho)"));
        surfaceScalarField rho_nei(fvc::interpolate(rho_, nei_, "reconstruct(rho)"));    

    surfaceVectorField rhoU_own(fvc::interpolate(rhoU_, own_, "reconstruct(U)"));
        surfaceVectorField rhoU_nei(fvc::interpolate(rhoU_, nei_, "reconstruct(U)"));

    surfaceVectorField U_own("U_own", rhoU_own/rho_own);
        surfaceVectorField U_nei("U_nei", rhoU_nei/rho_nei);

    volScalarField e("e", thermo_.h()-thermo_.p()/rho_); // -> stable
        //volScalarField e("e", rhoE_/(rho_)-0.5*magSqr(rhoU_/(rho_)));
        surfaceScalarField e_own(fvc::interpolate(e, own_, "reconstruct(T)"));
        surfaceScalarField e_nei(fvc::interpolate(e, nei_, "reconstruct(T)"));

    surfaceScalarField E_own("E_own", e_own + 0.5*magSqr(U_own));
        surfaceScalarField E_nei("E_nei", e_nei + 0.5*magSqr(U_nei));

    volScalarField rPsi("rPsi", 1.0/thermo_.psi());
        surfaceScalarField rPsi_own(fvc::interpolate(rPsi, own_, "reconstruct(T)"));
        surfaceScalarField rPsi_nei(fvc::interpolate(rPsi, nei_, "reconstruct(T)"));

    surfaceScalarField p_own("p_own", rho_own*rPsi_own);
        surfaceScalarField p_nei("p_nei", rho_nei*rPsi_nei);

    surfaceScalarField H_own(E_own+p_own/rho_own); 
    surfaceScalarField H_nei(E_nei+p_nei/rho_nei);         

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

    //covariant velocities
    surfaceVectorField normal("normal", mesh_.Sf()/mesh_.magSf());
    surfaceScalarField Uv_own=U_own & normal;
    surfaceScalarField Uv_nei=U_nei & normal;

    tmp <surfaceScalarField> S_own;
    tmp <surfaceScalarField> S_nei;

    if (directWaveSpeed)//direct wave speed estimates 
    {    
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
    
        surfaceScalarField UvTilde=UTilde() & normal;
        S_own=min(Uv_own-a_own, UvTilde-aTilde);
        S_nei=max(Uv_nei+a_nei, UvTilde+aTilde);
    }    

    tmp <surfaceScalarField> S_star
        (
               (p_nei-p_own+rho_own*Uv_own*(S_own()-Uv_own)-rho_nei*Uv_nei*(S_nei()-Uv_nei))/
        (rho_own*(S_own()-Uv_own)-rho_nei*(S_nei()-Uv_nei))
        );

        tmp <surfaceScalarField> pstar_own(p_own+rho_own*(S_own()-Uv_own)*(S_star()-Uv_own));
    tmp <surfaceScalarField> pstar_nei(p_nei+rho_nei*(S_nei()-Uv_nei)*(S_star()-Uv_nei)); 

    surfaceScalarField cf_own("cf_own",DST::pos(S_own()));                         // S_own>0 ? 1 : 0
    surfaceScalarField cf_m1("cf_m1",  DST::neg0(S_own())*DST::pos(S_star()));             // S_own<=0 and S_star>0 ? 1 : 0
    surfaceScalarField cf_m2("cf_m2",  DST::neg0(S_own())*DST::neg0(S_star())*DST::pos0(S_nei()));    // S_own<=0 and S_star<=0 and S_nei>=0 ? 1 : 0                          
        surfaceScalarField cf_nei("cf_nei",DST::neg0(S_own())*DST::neg0(S_star())*DST::neg(S_nei()));    // S_own<=0 and S_star<=0 and S_nei<0 and S_nei<=0 ? 1 : 0  

    surfaceScalarField rhoE_own("rhoE_own", rho_own*E_own);
    surfaceScalarField rhoE_nei("rhoE_nei", rho_nei*E_nei);    

    surfaceVectorField rhoUPhi_own ("rhoUPhi_own", rhoU_own*Uv_own+p_own*normal);
    surfaceVectorField rhoUPhi_nei ("rhoUPhi_nei", rhoU_nei*Uv_nei+p_nei*normal);

    surfaceScalarField rhoEPhi_own("rhoEPhi_own", (rhoE_own+p_own)*Uv_own);
    surfaceScalarField rhoEPhi_nei("rhoEPhi_nei", (rhoE_nei+p_nei)*Uv_nei);    

    phi_ = 
        (
        cf_own*rho_own*Uv_own +
              cf_m1*S_star()*rho_own*(S_own()-Uv_own)/(S_own()-S_star()) + 
        cf_m2*S_star()*rho_nei*(S_nei()-Uv_nei)/(S_nei()-S_star()) + 
              cf_nei*rho_nei*Uv_nei
    )*mesh_.magSf();

        phiUp_ = 
    (
        cf_own*rhoUPhi_own +
        cf_m1*(S_star()*(S_own()*rhoU_own-rhoUPhi_own)+S_own()*pstar_own()*normal)/(S_own()-S_star()) + 
        cf_m2*(S_star()*(S_nei()*rhoU_nei-rhoUPhi_nei)+S_nei()*pstar_nei()*normal)/(S_nei()-S_star()) + 
        cf_nei*rhoUPhi_nei
    )*mesh_.magSf(); 

    phiEp_ =
    (
        cf_own*rhoEPhi_own +
        cf_m1*S_star()*(S_own()*rhoE_own-rhoEPhi_own+S_own()*pstar_own())/(S_own()-S_star()) +
        cf_m2*S_star()*(S_nei()*rhoE_nei-rhoEPhi_nei+S_nei()*pstar_nei())/(S_nei()-S_star()) +
        cf_nei*rhoEPhi_nei
    )*mesh_.magSf();

    Up_ =     cf_own*U_own + 
        cf_m1*(S_own()-Uv_own)*U_own/(S_own()-S_star()) + 
        cf_m2*(S_nei()-Uv_nei)*U_nei/(S_nei()-S_star()) +
        cf_nei*U_nei; 

    calcFields();
        JacobianPtr_().computeJacobianDST();  
    
} // End of computeFluxDST() function

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
