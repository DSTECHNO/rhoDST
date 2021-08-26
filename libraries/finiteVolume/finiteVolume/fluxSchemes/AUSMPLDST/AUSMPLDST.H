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
    Foam::AUSMPLDST

Description
    Calculate fluxes using AUSM+ method

SourceFiles
    AUSMPLDST.C

References
-----------
    MS. Liu (1996). A sequel to AUSM: AUSM+. J Comp Physics, 
    129, 364-382.
    
\*-------------------------------------------------------------------------------------*/

#ifndef AUSMPLDST_H
#define AUSMPLDST_H

#include "baseFluxDST.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class AUSMPLDST Declaration
\*---------------------------------------------------------------------------*/

class AUSMPLDST 
: 
    public baseFluxDST
{
    protected:
    // Protected data

    // Private Data        
    private: 
        //- Coefficients
            const scalar alpha_ = 3.0/16.0;
            const scalar beta_ = 0.125;
    public:
        //- Runtime type information
            TypeName("AUSM+");

        //- Surface variable reference to owner cell
                surfaceScalarField own_;
        //- Surface variable reference to neighbour cell
        surfaceScalarField nei_;        

        //- public functions

        //- constructor
        AUSMPLDST
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
        );

        ~AUSMPLDST();

        virtual void computeFluxDST();        
                
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}// End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
