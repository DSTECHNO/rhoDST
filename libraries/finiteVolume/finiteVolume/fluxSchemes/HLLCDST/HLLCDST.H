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
    Foam::HLLCDST

Description
    Calculate fluxes using HLLCDST method

SourceFiles
    HLLCDST.C

References
-----------
    Eleuterio F. Toro, Michael Spruce, and William Speares. Restoration of the
     contact surface in the HLL-Riemann solver. Shock waves, 4(1):25–34, 1994.
    
\*---------------------------------------------------------------------------------------*/

#ifndef HLLCDST_H
#define HLLCDST_H

#include "baseFluxDST.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class HLLCDST Declaration
\*---------------------------------------------------------------------------*/

class HLLCDST 
: 
    public baseFluxDST
{
    // private data
    private: 
        bool directWaveSpeed; //Method for the calculation of wave speed 

    protected:
            // Protected data
            
    public:
        //- Runtime type information
            TypeName("HLLC");

        //- Surface variable reference to owner cell
            surfaceScalarField own_;
        //- Surface variable reference to neighbour cell
            surfaceScalarField nei_;

        //- public functions

        //- constructor
            HLLCDST
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

            ~HLLCDST();

            virtual void computeFluxDST();                
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
