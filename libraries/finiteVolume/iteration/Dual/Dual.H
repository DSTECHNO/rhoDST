/*-----------------------------------------------------------------------------------------*\
| =========                 |                                                               |
| \\      /  F ield         | rhoDST: Block-coupled solver for compressible flows           |
|  \\    /   O peration     |								    |
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
    rhoDST::Dual

Description
    Dual time iteration.    

SourceFiles
    Dual.C

\*---------------------------------------------------------------------------------------*/

#ifndef Dual_H
#define Dual_H

#include "iteration.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Dual Declaration
\*---------------------------------------------------------------------------*/

class Dual
:
    public iteration
{
    public:
        //- Runtime type information
            TypeName("Dual");

        //- Construct from components
            Dual
            (
                IOdictDST& dict,
                timeDST& runTime,
                basicThermo& thermo,
                volVectorField& U,
                compressible::turbulenceModel& turbulence,
                baseFluxDST& rhoDSTFlux
            );

        //- Destructor
            ~Dual();

        //- Member Functions
            virtual void iterate
            (
                volScalarField& rho,
                volVectorField& rhoU,
                volScalarField& rhoE
            );
    
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace rhoDST
#endif
