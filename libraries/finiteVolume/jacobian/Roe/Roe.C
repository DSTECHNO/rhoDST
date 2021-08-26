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

#include "JacobianDST.H"
#include "DST.H"
#include "Roe.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
defineNamedTemplateTypeNameAndDebug(Roe<tensor5>, 0);
addTemplatedToRunTimeSelectionTable(JacobianDST, Roe, tensor5, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Roe<Type>::Roe
(
    const fvMesh& mesh,
    const basicThermo& thermo,
    const compressible::turbulenceModel& turbulence,
    const volScalarField& rho,
    const volVectorField& rhoU,
    const volScalarField& rhoE
)
:
    ownerCoeffs_(mesh),
    neighbourCoeffs_(mesh),
    spectInternal_(mesh.nCells()),
    patchInternalCoeffs_(mesh.boundary().size()),
    patchNeighbourCoeffs_(mesh.boundary().size()),
    boundaryJacobian_(mesh.boundary().size()),
    spectPatch_(mesh.boundary().size()),
    mesh_(mesh),
    thermo_(thermo),
    turbulence_(turbulence),
    rho_(rho),
    rhoU_(rhoU),
    rhoE_(rhoE),
    nSf(mesh.Sf()/mesh.magSf()),
    Uf
    (
        IOobject
        (
            "Uf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("Uf",dimless,vector(0,0,0))
    ),
    gammaf
    (
        IOobject
        (
            "gammaf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("gammaf",dimless,scalar(0))
    ),
    phiJf
    (
        IOobject
        (
            "phiJf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("phiJf",dimless,scalar(0))
    ),
    a1f
    (
        IOobject
        (
            "a1f",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("a1f",dimless,scalar(0))
    )
{
    const IOdictionary& dict = mesh_.lookupObject<IOdictionary>("rhoDSTDict");
    addSpectralRadiusToPhysicalBoundaries = dict.lookupOrDefault("addSpectralRadiusToPhysicalBoundaries",false);
    forAll(mesh_.boundary(),patchI)
    {
       patchInternalCoeffs_[patchI] = Field<Type>(mesh_.Sf().boundaryField()[patchI].size(), Type::zero);
       patchNeighbourCoeffs_[patchI] = Field<Type>(mesh_.Sf().boundaryField()[patchI].size(), Type::zero);
       boundaryJacobian_[patchI] = Field<Type>(mesh_.Sf().boundaryField()[patchI].size(), Type::zero);
       spectPatch_[patchI] = Field<Type>(mesh_.Sf().boundaryField()[patchI].size(), Type::zero);
    } 
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Roe<Type>::~Roe()
{}

// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

template<class Type>
void Roe<Type>::computeJacobianDST()
{   

    // the notation and formulations are taken from Blazek(2001) A.7 pp. 419-420 
    volVectorField U_(rhoU_/rho_);
    volScalarField gamma(thermo_.Cp()/thermo_.Cv());
    volScalarField phiJ(0.5*(gamma-1)*magSqr(U_));
    volScalarField a1(gamma*rhoE_/rho_-phiJ);
    
    // ownerCoeffs: a surfaceTypeField that stores face owner's jabobians.
    calcJ(ownerCoeffs_, true, U_, gamma, phiJ, a1, mesh_.owner());
    
    // neighbourCoeffs: a surfaceTypeField that stores face owner's jabobians.
    calcJ(neighbourCoeffs_, false, U_, gamma, phiJ, a1, mesh_.neighbour());

    ownerCoeffs_ *= mesh_.weights();
    neighbourCoeffs_ *= (1-mesh_.weights());
    
    forAll(mesh_.boundary(),patchI)
    {
        const scalarField magSf = mesh_.magSf().boundaryField()[patchI];
        if(mesh_.boundary()[patchI].coupled())
        {
            patchInternalCoeffs_[patchI] = ownerCoeffs_.boundaryField()[patchI]*magSf;
            patchNeighbourCoeffs_[patchI] = neighbourCoeffs_.boundaryField()[patchI]*magSf;
        }        
    }

    surfaceScalarField amaxSf
    (
        mag(fvc::interpolate(U_)&mesh_.Sf()) +
        mesh_.magSf()*fvc::interpolate(sqrt(thermo_.Cp()/thermo_.Cv()/thermo_.psi()))
    );

    surfaceScalarField muEffSf("muEffSf", fvc::interpolate(turbulence_.muEff()));

    surfaceScalarField alphaEff("alphaEff", fvc::interpolate(turbulence_.alphaEff()));

    surfaceScalarField rhoSf("rhoSf", fvc::interpolate(thermo_.rho()));
    
    Type I = DST::tensorI<Type>();
    
    const surfaceScalarField& deltaC =  mesh_.deltaCoeffs();

    const surfaceScalarField& magSf =  mesh_.magSf();

    surfaceScalarField lambda = 0.5*(amaxSf + magSf*(muEffSf+alphaEff)/rhoSf*deltaC);

    spectInternal_ = I*lambda.internalField();
    
    forAll(mesh_.boundary(),patchI)
    {
        // Spectral radius is added to the coupled patches only.
        if(mesh_.boundary()[patchI].coupled())
        {
            spectPatch_[patchI] = I*(lambda.boundaryField()[patchI]);
        }
        // TODO: cyclic patches must be considered
    }  

    computeBoundaryJacobians();  
}

template<class Type>
template<class fType>
void Roe<Type>::updateCoeffs
(
    GeometricField<fType, fvsPatchField, surfaceMesh>& fJ,
    const GeometricField<fType, fvPatchField, volMesh> vJ,
    const unallocLabelList& addr,
    bool owner
)
{
    forAll(fJ,faceI)
    {
        fJ[faceI] = vJ[addr[faceI]];
    }
    
    forAll(fJ.boundaryField(), patchI)
    {
        if(fJ.boundaryField()[patchI].coupled())
        {
            if(owner)
            {
                fJ.boundaryField()[patchI] = vJ.boundaryField()[patchI].patchInternalField();
            }
            else
            {
                fJ.boundaryField()[patchI] = vJ.boundaryField()[patchI].patchNeighbourField();
            }
        }
    }
}

template<class Type>
void Roe<Type>::calcJ
(
    jacobian<Type>& J,
    bool owner,
    const volVectorField& U,
    const volScalarField& gamma,
    const volScalarField& phiJ,
    const volScalarField& a1,
    const unallocLabelList& addr
)
{
    updateCoeffs(Uf, U, addr, owner);
    updateCoeffs(gammaf, gamma, addr, owner);
    updateCoeffs(phiJf, phiJ, addr, owner);
    updateCoeffs(a1f, a1, addr, owner);
    

    J.insertScalarEqn(nSf ,0,1);
    J.insertVectorEqn((phiJf*nSf - ((Uf*Uf)&nSf)),1,0);
    J.insertVectorEqn(Uf*nSf + (Uf&nSf)*tensor::I -(gammaf-1)*(nSf*Uf),1,1);
    J.insertVectorEqn((gammaf-1)*nSf,1,4);
    J.insertScalarEqn((Uf&nSf)*(phiJf-a1f),4,0);
    J.insertScalarEqn(a1f*nSf - (Uf&nSf)*(gammaf -1)*Uf,4,1);
    J.insertScalarEqn(gammaf*(Uf&nSf),4,4);
}

template<class Type>
void Roe<Type>::computeBoundaryJacobians()
{ 
    const fvMesh& mesh(mesh_);

    const volScalarField& rho(rho_);
    const volVectorField& rhoU(rhoU_);
    const volScalarField& rhoE(rhoE_);
    const volScalarField& p(thermo_.p());
    const volScalarField& T(thermo_.T());
    const volScalarField& h(thermo_.h()); 

    volVectorField U_(rhoU_/rho_);
    
    /*---------------------------------------------------------------------
        
    \Delta R_b =\left ( \frac{\partial R}{\partial P} \right )_b
                \left ( \frac{\partial P_b}{\partial P_i}  \right )
                \left ( \frac{\partial P}{\partial W}  \right )_i\Delta W_i
    
    Here "b" and "i" subscripts denote the variable at boundaries 
    and at neighbouring cell centers respectively:
    
    W_i = {rho, rhoU, rhoE}
    W_b = {rho_b, rhoU_b, rhoE_b}
    P_i = {p, U, T}
    P_b = {p_b, U_b, T_b}
    
    
    ---------------------------------------------------------------------*/
    
    
    forAll(mesh.boundary(),patchI)
    {
        if(!mesh.boundary()[patchI].coupled())
        {
            const fvPatch& patch = mesh.boundary()[patchI];
            const vectorField& pSf = patch.Sf();

            //calculate (dR/dX) on the physical boundary
            Field<Type> dRdXb(pSf.size(),Type::zero);

            vectorField UB = rhoU.boundaryField()[patchI]/rho.boundaryField()[patchI];

            const scalarField& pB = p.boundaryField()[patchI];
            const scalarField& TB = T.boundaryField()[patchI];
            const scalarField& rhoB = rho.boundaryField()[patchI];
            scalarField CvB = thermo_.Cv()().boundaryField()[patchI];

            scalarField EB = rhoE.boundaryField()[patchI]/rho.boundaryField()[patchI];
            scalarField rRT = rhoB/pB;

            DST::insertScalarEqn(dRdXb, (UB & pSf)*rRT, 0, 0);
            DST::insertScalarEqn(dRdXb, rhoB*pSf, 0, 1);
            DST::insertScalarEqn(dRdXb, -rhoB/TB*(UB & pSf), 0, 4);
            DST::insertVectorEqn(dRdXb, rRT*(UB & pSf)*UB+pSf,1,0);
            DST::insertVectorEqn(dRdXb, rhoB*UB*pSf + tensor::I*rhoB*(UB & pSf),1,1);
            DST::insertVectorEqn(dRdXb, -rhoB/TB*(UB & pSf)*UB,1,4);
            DST::insertScalarEqn(dRdXb, (EB*rRT+1)*(UB & pSf), 4, 0);
            DST::insertScalarEqn(dRdXb, rhoB*(UB & pSf)*UB+pSf*(rhoB*EB+pB), 4, 1);
            DST::insertScalarEqn(dRdXb, rhoB*CvB*(UB & pSf), 4, 4);


            //calculate (dX/dW) at a patch internal field
            vectorField UI = U_.boundaryField()[patchI].patchInternalField();
            scalarField rhoEI = rhoE_.boundaryField()[patchI].patchInternalField();

            scalarField p = thermo_.p().boundaryField()[patchI].patchInternalField();
            scalarField rRho = scalar(1)/rho_.boundaryField()[patchI].patchInternalField();
            scalarField Cv = thermo_.Cv()().boundaryField()[patchI].patchInternalField();

            scalarField rRhoCv = rRho/Cv;
            scalarField e = rhoEI*rRho - 0.5*magSqr(UI);
            scalarField RT = p*rRho;
            scalarField RCv = RT/e;            

            const scalarField& weights = mesh.weights().boundaryField()[patchI];

            Field<Type> dXidWi(weights.size(),Type::zero);
            
            DST::insertScalarEqn(dXidWi, RT-RCv*(e-0.5*magSqr(UI)), 0, 0);
            DST::insertScalarEqn(dXidWi, -RCv*UI, 0, 1);
            DST::insertScalarEqn(dXidWi, RCv, 0, 4);
            DST::insertVectorEqn(dXidWi, -rRho*UI,1,0);
            DST::insertVectorEqn(dXidWi, tensor::I*rRho,1,1);
            DST::insertScalarEqn(dXidWi, rRhoCv*(0.5*magSqr(UI)-e), 4, 0);
            DST::insertScalarEqn(dXidWi, -rRhoCv*UI, 4, 1);
            DST::insertScalarEqn(dXidWi, rRhoCv, 4, 4);

            //calculate (dX_b)/(dX_i) on the boundary                   
            const vectorField Uf = (mesh.lookupObject<volVectorField>("U")).boundaryField()[patchI].valueInternalCoeffs(weights);
            const scalarField pf = thermo_.p().boundaryField()[patchI].valueInternalCoeffs(weights);
            const scalarField Tf = thermo_.T().boundaryField()[patchI].valueInternalCoeffs(weights);
            
            Field<Type> dXbdXi(weights.size(),Type::zero);

            DST::insertScalarEqn(dXbdXi, pf, 0, 0);
            DST::insertScalarEqn(dXbdXi, Uf.component(0), 1, 1);
            DST::insertScalarEqn(dXbdXi, Uf.component(1), 2, 2);
            DST::insertScalarEqn(dXbdXi, Uf.component(2), 3, 3);
            DST::insertScalarEqn(dXbdXi, Tf, 4, 4);

            patchInternalCoeffs_[patchI] = ((dRdXb & dXbdXi) & dXidWi);
        }
    }
}

template<class Type>
jacobian<Type>& Roe<Type>::ownerCoeffs()
{
    return ownerCoeffs_;
}
template<class Type>
jacobian<Type>& Roe<Type>::neighbourCoeffs()
{
    return neighbourCoeffs_;
}
template<class Type>
Field<Field<Type>>& Roe<Type>::patchInternalCoeffs()
{
    return patchInternalCoeffs_;
}
template<class Type>
Field<Field<Type>>& Roe<Type>::patchNeighbourCoeffs()
{
    return patchNeighbourCoeffs_;
}

template<class Type>
Field<Field<Type>>& Roe<Type>::boundaryJacobian()
{
    return boundaryJacobian_;
}

template<class Type>
Field<Type>& Roe<Type>::spectInternal()
{
    return spectInternal_;
}
template<class Type>
Field<Field<Type>>& Roe<Type>::spectPatch()
{
    return spectPatch_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam






