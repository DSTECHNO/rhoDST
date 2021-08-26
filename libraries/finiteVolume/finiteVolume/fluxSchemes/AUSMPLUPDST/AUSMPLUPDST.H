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
    Foam::AUSMPLUPDST

Description
    Calculate fluxes using AUSM+up method

SourceFiles
    AUSMPLUPDST.C

References
-----------
    MS. Liu (2006). A sequel to AUSM, Part II: AUSM+-up for all speeds. 
    J Comp Physics, 214, 137-170
    
\*----------------------------------------------------------------------------------------*/

#ifndef AUSMPLUPDST_H
#define AUSMPLUPDST_H

#include "baseFluxDST.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class AUSMPLUPDST Declaration
\*---------------------------------------------------------------------------*/


class AUSMPLUPDST 
: 
    public baseFluxDST
{
    tmp <surfaceScalarField> M0;
    tmp <surfaceScalarField> Mmean2;
    tmp <surfaceScalarField> fa;
    tmp <surfaceScalarField> alpha;    

    protected:
            // Protected data
            
    // Private Data        
    private: 
        //- Empirical coefficients suggested by Liou        
            const scalar beta_ = 0.125; 
        const scalar Kp_ = 0.25; 
        const scalar Ku_ = 0.75;
        const scalar sigma_=1.0;

        //- low Mach number dissipation in AUSM+up
        scalar Mainf;
        bool lowMach;

    public:
        //- Runtime type information
            TypeName("AUSM+up");

        //- Surface variable reference to owner cell
                surfaceScalarField own_;
        //- Surface variable reference to neighbour cell
        surfaceScalarField nei_;

        //- public functions

        //- constructor
        AUSMPLUPDST
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

        ~AUSMPLUPDST();

        void ReadLowMachDiss();

        virtual void computeFluxDST();        
                
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}// End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
