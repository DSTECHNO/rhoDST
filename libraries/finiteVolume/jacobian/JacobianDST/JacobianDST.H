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
    rhoDST::JacobianDST

Description
     Calculates Jacobians using analytical descriptions given in the following reference. 

SourceFiles
    JacobianDST.C

References
-----------
    J Blazek. Computational Fluid Dynamics: Principles and Applications.
    
\*---------------------------------------------------------------------------------------*/

#ifndef JacobianDST_H
#define JacobianDST_H

#include "IOdictionary.H"
#include "IOdictDST.H"
#include "volFields.H"
#include "fvCFD.H"
#include "Vector5typedefs.H"
#include "basicPsiThermo.H"
#include "turbulenceModel.H"
#include "jacobian.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class JacobianDST Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class JacobianDST
{
    private:
        JacobianDST(const JacobianDST&);
            void operator=(const JacobianDST&);
    
    public:
            //- Runtime type information
            TypeName("JacobianDST");

            //- Declare run-time constructor selection table
            declareRunTimeSelectionTable
            (
                autoPtr,
                JacobianDST,
                dictionary,
                (
                    const fvMesh& mesh,
                    const basicThermo& thermo,
                    const compressible::turbulenceModel& turbulence,
                    const volScalarField& rho,
                    const volVectorField& rhoU,
                    const volScalarField& rhoE
                ),
                (mesh, thermo, turbulence, rho, rhoU, rhoE)
            );
     
            //- Constructors
            JacobianDST();

        //- Destructor
            virtual ~JacobianDST();

            //- Selector
            static autoPtr<JacobianDST<Type>> New
            (
                const IOdictionary& dict,
                const fvMesh& mesh,
                const basicThermo& thermo,
                const compressible::turbulenceModel& turbulence,
                const volScalarField& rho,
                const volVectorField& rhoU,
                const volScalarField& rhoE
            );    
                
        //- Member Functions
            virtual void computeJacobianDST()=0;
            virtual jacobian<Type>& ownerCoeffs()=0;
            virtual jacobian<Type>& neighbourCoeffs()=0;
            virtual Field<Field<Type>>& patchInternalCoeffs()=0;
            virtual Field<Field<Type>>& patchNeighbourCoeffs()=0;
            virtual Field<Field<Type>>& boundaryJacobian()=0;
            virtual Field<Type>& spectInternal()=0;
            virtual Field<Field<Type>>& spectPatch()=0;
    
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam


#ifdef NoRepository
#   include "JacobianDST.C"
#endif

#endif
