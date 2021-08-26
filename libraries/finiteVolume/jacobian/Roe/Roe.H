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
                                      

Class
    rhoDST::Roe

Description
    Roe flux splitting using mesh weights.
    R = (nF_{i}+(1-n)F_{j}) - \frac{1}{2}\left | \lambda \right |(U_{j}-U_{i})

SourceFiles
    Roe.C

\*---------------------------------------------------------------------------------------*/

#ifndef Roe_H
#define Roe_H

#include "IOdictionary.H"
#include "volFields.H"
#include "fvCFD.H"
#include "basicPsiThermo.H"
#include "JacobianDST.H"


namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class JacobianDST Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Roe 
: 
public JacobianDST<Type>
{
    protected:
        jacobian<Type> ownerCoeffs_;
        jacobian<Type> neighbourCoeffs_;
        Field<Type> spectInternal_;
        Field<Field<Type>> patchInternalCoeffs_;
        Field<Field<Type>> patchNeighbourCoeffs_;
        Field<Field<Type>> boundaryJacobian_;
        Field<Field<Type>> spectPatch_;

        const fvMesh& mesh_;
        const basicThermo& thermo_;
        const compressible::turbulenceModel& turbulence_;
        const volScalarField& rho_;
        const volVectorField& rhoU_;
        const volScalarField& rhoE_;

        surfaceVectorField nSf;
        surfaceVectorField Uf;
        surfaceScalarField gammaf;
        surfaceScalarField phiJf;
        surfaceScalarField a1f;
        
        bool addSpectralRadiusToPhysicalBoundaries;

    public:
        //- Runtime type information
        TypeName("Roe");

        //- Construct from components
        Roe
        (
            const fvMesh& mesh,
            const basicThermo& thermo,
            const compressible::turbulenceModel& turbulence,
            const volScalarField& rho,
            const volVectorField& rhoU,
            const volScalarField& rhoE
        );

        //- Destructor
        virtual ~Roe();

        //- Member Functions
        virtual void computeJacobianDST();    
    
        void computeBoundaryJacobians();
    
        template<class fType>
        void updateCoeffs
        (
            GeometricField<fType, fvsPatchField, surfaceMesh>& fJ,
            const GeometricField<fType, fvPatchField, volMesh> vJ,
            const unallocLabelList& addr,
            bool owner
        );
    
        void calcJ
        (
            jacobian<Type>& J,
            bool owner,
            const volVectorField& U,
            const volScalarField& gamma,
            const volScalarField& phiJ,
            const volScalarField& a1,
            const unallocLabelList& addr
        );
    
        virtual jacobian<Type>& ownerCoeffs();
        virtual jacobian<Type>& neighbourCoeffs();
        virtual Field<Field<Type>>& patchInternalCoeffs();
        virtual Field<Field<Type>>& patchNeighbourCoeffs();
        virtual Field<Field<Type>>& boundaryJacobian();
        virtual Field<Type>& spectInternal();
        virtual Field<Field<Type>>& spectPatch();
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace rhoDST
#endif
