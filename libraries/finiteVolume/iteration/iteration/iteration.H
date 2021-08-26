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
    rhoDST::iteration

Description
    Iteration base class for rhoDST solver.

SourceFiles
    iteration.C

\*---------------------------------------------------------------------------------------*/

#ifndef iteration_H
#define iteration_H

#include "IOdictionary.H"
#include "IOdictDST.H"
#include "volFields.H"
#include "fvCFD.H"
#include "Vector5typedefs.H"
#include "fvBlockMatrix.H"
#include "basicPsiThermo.H"
#include "turbulenceModel.H"
#include "fvBlockMatrixDST.H"
#include "baseFluxDST.H"
#include "JacobianDST.H"
#include "fvmDST.H"
#include "DST.H"
#include "bound.H"
#include "timeDST.H"
#include "geometricFieldDST.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class iteration Declaration
\*---------------------------------------------------------------------------*/

class iteration 
{
private:
    //- Construct as copy (not implemented)
    iteration(const iteration&);

    //- Disallow default bitwise assignment
        void operator=(const iteration&);        

    protected:        
            timeDST& runTime_;
            basicThermo& thermo_;
            volVectorField& U_;
            const fvMesh& mesh_;
            compressible::turbulenceModel& turbulence_;
            baseFluxDST& rhoDSTFlux_;
            volVector5FieldDST W;
    
    public:
            //- Runtime type information
            TypeName("iteration");

            //- Declare run-time constructor selection table
            declareRunTimeSelectionTable
            (
                autoPtr,
                iteration,
                dictionary,
                (
                    const IOdictDST& dict,
                    timeDST& runTime,
                    basicThermo& thermo,
                    volVectorField& U,
                    compressible::turbulenceModel& turbulence,
                    baseFluxDST& rhoDSTFlux
                ),
                (dict, runTime, thermo, U, turbulence, rhoDSTFlux)
               );

            //- Construct from components
            iteration
            (
                const IOdictDST& dict,
                timeDST& runTime,
                basicThermo& thermo,
                volVectorField& U,
                compressible::turbulenceModel& turbulence,
                baseFluxDST& rhoDSTFlux
            );

            //- Destructor
            virtual ~iteration();

            //- Selector
            static autoPtr<iteration> New
            (
                    const IOdictDST& dict,
                    timeDST& runTime,
                    basicThermo& thermo,
                    volVectorField& U,
                    compressible::turbulenceModel& turbulence,
                    baseFluxDST& rhoDSTFlux
            );        

            //- Member Functions
            volVector5FieldDST& getW();

            virtual void iterate
            (
                    volScalarField& rho,
                    volVectorField& rhoU,
                    volScalarField& rhoE
            )=0;
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace rhoDST

#endif
