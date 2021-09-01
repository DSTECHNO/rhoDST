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
    Foam::HLLCPDST

Description
    Calculate fluxes using HLLCPDST method

SourceFiles
    HLLCPDST.C

References
-----------
    Xie W., Lai J. and Li.H. An accurate and roboust HLLC-type Riemann solver 
    for the    compressible Euler system at various Mach numbers. 
    Int. J. for Num. Fluids, 430-463, 2008.
    
\*---------------------------------------------------------------------------------------*/

#ifndef HLLCPDST_H
#define HLLCPDST_H

#include "baseFluxDST.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class HLLCPDST Declaration
\*---------------------------------------------------------------------------*/

class HLLCPDST 
: 
    public baseFluxDST
{
    // private data
    private: 
        bool directWaveSpeed; //Method for the calculation of wave speed
        tmp<surfaceScalarField> f; 

    protected:
            // Protected data
            
    public:
        //- Runtime type information
            TypeName("HLLCP");

        //- Surface variable reference to owner cell
            surfaceScalarField own_;
        //- Surface variable reference to neighbour cell
            surfaceScalarField nei_;

        //- public functions

            HLLCPDST
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
            );

            ~HLLCPDST();

            virtual void computeFluxDST();  

        //- calculate f falues considering neighboring cells
            virtual void computef(const volScalarField& p);              
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
