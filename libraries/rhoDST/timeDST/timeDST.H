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
    rhoDST::timeDST

Description
    Time class for rhoDST

SourceFiles
    timeDST.C

   
\*---------------------------------------------------------------------------*/

#ifndef timeDST_H
#define timeDST_H

#include "IOdictionary.H"
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

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class timeDST Declaration
\*---------------------------------------------------------------------------*/
class timeDST
:
    public Time
{
    private:

    protected:
            basicThermo* thermoPtr;
            compressible::turbulenceModel* turbulencePtr;
            volVectorField* UPtr;
            const IOdictDST* dictPtr;
            autoPtr<volScalarField> rhoOldPtr;
            autoPtr<volScalarField> eOldPtr;
            scalar maxCo;
            scalar minCo;
            scalar courant;
            autoPtr<volScalarField> courantPtr;
            autoPtr<volScalarField> rDeltaTDSTPtr;
            scalar residual;
            scalar residualRatio;

            // blazek courant
            autoPtr<volVectorField> dS;
            bool neighbouring;

            int innerIter;
            int maxInnerIter;        
                
    public:
            //- Constructors
            timeDST
            (
                    const argList& args
            );

            //- Destructor
            ~timeDST();

            //- Member Functions
            void initialize
            (
                basicThermo& thermo,
                const IOdictDST& dict,
                compressible::turbulenceModel& turbulence,
                volVectorField& U
            );

            bool correctTurbulence();
            bool steadyState();
            scalar coNum();
            void setTransientDeltaT();
            void setSteadyStateDeltaT();
            void setCourantDST(const BlockSolverPerformance<vector5>* solverPerf = nullptr);
            bool localTimeStepping();
            const scalarField& rDeltaTDST() const;
            int nIter() const;
            bool innerLoop();
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace rhoDST
#endif
