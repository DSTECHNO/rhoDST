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
    rhoDST::KurganovDST

Description
    Calculate fluxes using Kurganov method

SourceFiles
    KurganovDST.C

References:
----------- 
    A. Kurganov and E. Tadmor (2000). New high-resolution central schemes for nonlinear 
    conservation and convection-diffusion equations

    C.J. Greenshields, H.G. Weller, L. Gasparini and J.M. Reese (2009). Implementation of 
    semi-discrete, non-staggered central scehemes in a colocated, polyhedral, finite volume
    framework, for high-speed viscous flows. 
    
\*--------------------------------------------------------------------------------------------*/                  
    
#ifndef KurganovDST_H
#define KurganovDST_H

#include "baseFluxDST.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class KurganovDST Declaration
\*---------------------------------------------------------------------------*/

class KurganovDST 
: 
    public baseFluxDST
{

    protected:

    public:
            //- Runtime type information
            TypeName("Kurganov");

           //- Dimensionless scalar variable for zero
            dimensionedScalar v_zero;

            //- constructor
            KurganovDST
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

            //- Destructor
            ~KurganovDST();
        
            //- Member functions
            virtual void computeFluxDST();                
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
}// End namespace Foam
#endif
