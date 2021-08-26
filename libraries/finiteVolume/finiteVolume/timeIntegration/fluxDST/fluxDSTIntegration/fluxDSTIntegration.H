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
    Foam::fluxDSTIntegration

Description
    Base class for the time integration of governing equations

SourceFiles
    fluxDSTIntegration.C

\*---------------------------------------------------------------------------------------*/

#ifndef fluxDSTIntegration_H
#define fluxDSTIntegration_H

#include "IOdictionary.H"
#include "IOdictDST.H"
#include "volFields.H"
#include "fvCFD.H"
#include "Vector5typedefs.H"
#include "fvBlockMatrixDST.H"
#include "JacobianDST.H"
#include "runTimeSelectionTables.H"
#include "geometricFieldDST.H"

namespace Foam
{

template<class Type>
class geometricFieldDST;

namespace fvDST
{


/*---------------------------------------------------------------------------*\
                           Class fluxDST Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fluxDSTIntegration : public refCount // Base class
{
protected:
    //- Disallow default bitwise copy construct 
        fluxDSTIntegration(const fluxDSTIntegration&);

    //- Disallow default bitwise assignment
        void operator=(const fluxDSTIntegration&);


public:
    //- Runtime type information
        TypeName("fluxDSTIntegration");

    //- Declare run-time constructor selection table
        declareRunTimeSelectionTable
        (
            tmp,
            fluxDSTIntegration,
            dictionary,
            (),
            ()
        );
        
    //- Construct from components
        fluxDSTIntegration();

    //- destructor
        virtual ~fluxDSTIntegration();

    //- member functions
        static tmp<fluxDSTIntegration<Type>> New
        (
            const IOdictDST& dict
        );    
        
        virtual tmp<fvBlockMatrixDST<Type>> integrate
        (
            baseFluxDST& rhoDSTFlux,
            geometricFieldDST<Type>& W
        ) const=0; 
            
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
}

}

#ifdef NoRepository
#   include "fluxDSTIntegration.C"
#endif


#endif
