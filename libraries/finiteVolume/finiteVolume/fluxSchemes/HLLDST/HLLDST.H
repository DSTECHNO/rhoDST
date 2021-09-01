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
    Foam::HLLDST

Description
    Calculate fluxes using HLL method

SourceFiles
    HLLDST.C

References
-----------
    Eleuterio F. Toro, Michael Spruce, and William Speares. Restoration of the
     contact surface in the HLL-Riemann solver. Shock waves, 4(1):25–34, 1994.
    
\*---------------------------------------------------------------------------------------*/

#ifndef HLLDST_H
#define HLLDST_H

#include "baseFluxDST.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class HLLDST Declaration
\*---------------------------------------------------------------------------*/

class HLLDST 
: 
    public baseFluxDST
{

    protected:

    public:
        //- Runtime type information
            TypeName("HLL");

        //- constructor
        HLLDST
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

        ~HLLDST();

        virtual void computeFluxDST();                
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam
#endif
