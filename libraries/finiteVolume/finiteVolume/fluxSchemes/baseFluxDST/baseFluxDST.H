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
    Foam::baseFluxDST

Description
    Base class for flux schemes to interpolate fields 

SourceFiles
    baseFluxDST.C
    
\*-----------------------------------------------------------------------------------------*/

#ifndef baseFluxDST_H
#define baseFluxDST_H

#include "IOdictionary.H"
#include "IOdictDST.H"
#include "volFields.H"
#include "fvCFD.H"
#include "Vector5typedefs.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "basicPsiThermo.H"
#include "turbulenceModel.H"
#include "fvBlockMatrix.H"
#include "DST.H"
#include "JacobianDST.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class baseFluxDST Declaration
\*---------------------------------------------------------------------------*/

class baseFluxDST  // Base class
{
private:
    //- Construct as copy (not implemented)
        baseFluxDST(const baseFluxDST&);

    //- Disallow default bitwise assignment
        void operator=(const baseFluxDST&);        

protected:

    //- Jacobian pointer
        autoPtr<JacobianDST<tensor5>> JacobianPtr_;
    
    //- Whether viscous fluxes are added. Used in calcFields().
        bool inviscid_;
    
    //- References
        const fvMesh& mesh_;
        const basicThermo& thermo_;
        const volScalarField& rho_;
        const volVectorField& U_;
        const volVectorField& rhoU_;
        const volScalarField& rhoE_;
        const compressible::turbulenceModel& turbulence_;
    
    //- Density, Velocity and Energy Fluxes
        surfaceScalarField& phi_;
        surfaceVectorField phiUp_;
        surfaceScalarField phiEp_;
        surfaceVectorField Up_;
        
    //- Residual 
        volVector5Field R_;


    //- Surface variable reference to owner cell
        surfaceScalarField own_;
    
    //- Surface variable reference to neighbour cell
        surfaceScalarField nei_;

public:
    //- Runtime type information
        TypeName("baseFluxDST");

    //- Declare run-time constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            baseFluxDST,
            dictionary,
            (
                const IOdictionary& dict,
                const fvMesh& mesh,
                basicPsiThermo& thermo,
                const volScalarField& rho,
                volVectorField& U,
                const volVectorField& rhoU,
                const volScalarField& rhoE,
                compressible::turbulenceModel& turbulence,
                surfaceScalarField& phi
            ),
            (dict, mesh, thermo, rho, U, rhoU, rhoE, turbulence, phi)
        );

    //- Construct from components
        baseFluxDST
        (
            const word& type,
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
        virtual ~baseFluxDST();

    //- Selector
        static autoPtr<baseFluxDST> New
        (
            const IOdictionary& dict,
            const fvMesh& mesh,
            basicPsiThermo& thermo,
            const volScalarField& rho,
            volVectorField& U,
            const volVectorField& rhoU,
            const volScalarField& rhoE,
            compressible::turbulenceModel& turbulence,
            surfaceScalarField& phi
        );        

    //- Member functions
        virtual void computeFluxDST()=0;
    
        surfaceScalarField& phi();

        surfaceVectorField& phiUp();

        surfaceScalarField& phiEp();

        surfaceVectorField& Up();

        bool& inviscid();

        void calcFields();

        const fvMesh& mesh();

        JacobianDST<tensor5>& Jacobian();

        vector5Field& R();

        vector5Field& RoldTime();
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

#endif
