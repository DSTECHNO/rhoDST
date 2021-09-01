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
    rhoDST::geometricFieldDST

SourceFiles
    geometricFieldDST.C
  
\*---------------------------------------------------------------------------------------*/

#ifndef geometricFieldDST_H
#define geometricFieldDST_H

#include "IOdictionary.H"
#include "volFields.H"
#include "fvCFD.H"
#include "Vector5typedefs.H"
#include "timeDST.H"

namespace Foam
{

class timeDST;
/*---------------------------------------------------------------------------*\
                           Class geometricFieldDST Declaration
\*---------------------------------------------------------------------------*/
template<class Type>
class geometricFieldDST
:
    public GeometricField<Type, Foam::fvPatchField, Foam::volMesh>
{
    private:

    protected:
            GeometricField<Type, Foam::fvPatchField, Foam::volMesh>* deltaDST_;
            volScalarField* rho_;
            volVectorField* rhoU_;
            volScalarField* rhoE_;
            timeDST& runTimeDST_;
            IOdictDST& dictDST_;
    
    public:
            //- Constructors
                geometricFieldDST
                (
                    const IOobject& io,
                    const fvMesh& mesh,
                    const dimensioned<Type>& dt,
                    timeDST& runTimeDST,
                    IOdictDST& dictDST
                );

            //- Destructor
                ~geometricFieldDST();


            //- Member Functions
                timeDST& runTime();
                IOdictDST& dictDST();
                GeometricField<Type, Foam::fvPatchField, Foam::volMesh>& delta();

                template<class fType>
                void insert(GeometricField<fType, Foam::fvPatchField, Foam::volMesh>& f, label cmpt);

                template<class fType>
                void retrieve(GeometricField<fType, Foam::fvPatchField, Foam::volMesh>& f, label cmpt);
                void set
                (
                    volScalarField& rho,
                    volVectorField& rhoU,
                    volScalarField& rhoE
                );

                volScalarField& rho();
                volVectorField& rhoU();
                volScalarField& rhoE();

            //- Operators
                void operator+=
                (
                    const geometricFieldDST<Type>& gf
                );

                void operator+=
                (
                    const GeometricField<Type, Foam::fvPatchField, Foam::volMesh>& gf
                );

                geometricFieldDST<Type>& operator++(int);
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef geometricFieldDST<vector5> volVector5FieldDST;

}


#ifdef NoRepository
#   include "geometricFieldDST.C"
#endif

#endif
