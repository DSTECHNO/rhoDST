/*-----------------------------------------------------------------------------------------*\
| =========                 |                                                               |
| \\      /  F ield         | rhoDST: Block-coupled solver for compressible flows           |
|  \\    /   O peration     |                                                               |
|   \\  /    A nd           | Copyright (C) 2020 Design and Simulation Tech. Inc. (DSTECH)  |
|    \\/     M anipulation  | Website:  http://www.dstechno.net/                            |
\*-----------------------------------------------------------------------------------------*/
/*
License
    This file is derivative work of rhoDST.

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

\*---------------------------------------------------------------------------*/

#include "geometricFieldDST.H"
namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
geometricFieldDST<Type>::geometricFieldDST
(
    const IOobject& io,
    const fvMesh& mesh,
    const dimensioned<Type>& dt,
    timeDST& runTimeDST,
    IOdictDST& dictDST
)
:
    GeometricField<Type, Foam::fvPatchField, Foam::volMesh>
    (
        io,
        mesh,
        dt
    ),
    deltaDST_(NULL),
    runTimeDST_(runTimeDST),
    dictDST_(dictDST)
{} 

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
geometricFieldDST<Type>::~geometricFieldDST()
{
    if(deltaDST_)
        delete deltaDST_;
}


// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

template<class Type>
timeDST& geometricFieldDST<Type>::runTime()
{
    return runTimeDST_;
}

template<class Type>
GeometricField<Type, Foam::fvPatchField, Foam::volMesh>& geometricFieldDST<Type>::delta()
{
    if(!deltaDST_)
    {
        deltaDST_ = 
        new GeometricField<Type, Foam::fvPatchField, Foam::volMesh>
        (
            IOobject
            (
                "deltaW",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensioned<Type>("deltaDST_",dimless,Type::zero)
        );
    }
    return *deltaDST_;
}

template<class Type>
IOdictDST& geometricFieldDST<Type>::dictDST()
{
    return dictDST_;
}



template<class Type>
volScalarField& geometricFieldDST<Type>::rho()
{
    return *rho_;
}

template<class Type>
volVectorField& geometricFieldDST<Type>::rhoU()
{
    return *rhoU_;
}

template<class Type>
volScalarField& geometricFieldDST<Type>::rhoE()
{
    return *rhoE_;
}


template<class Type>
void geometricFieldDST<Type>::set
(
    volScalarField& rho,
    volVectorField& rhoU,
    volScalarField& rhoE
)
{
    rho_ = &rho;
    rhoU_ = &rhoU;
    rhoE_ = &rhoE;
    
    this->insert(rho, 0);
    this->insert(rhoU, 1);
    this->insert(rhoE, 4);
}



template<class Type>
template<class fType>
void geometricFieldDST<Type>::insert
(
    GeometricField<fType, Foam::fvPatchField, Foam::volMesh>& f, 
    label cmpt
)
{
    int nComponents = pTraits<fType>::nComponents;
    for(int i=0; i<nComponents; i++)
        this->replace(i+cmpt,f.component(i));
}

template<class Type>
template<class fType>
void geometricFieldDST<Type>::retrieve
(
    GeometricField<fType, Foam::fvPatchField, Foam::volMesh>& f, 
    label cmpt
)
{
    int nComponents = pTraits<fType>::nComponents;
    for(int i=0; i<nComponents; i++)
        f.replace(i,this->component(i+cmpt));
}
// * * * * * * * * * * * * * * * * Operators  * * * * * * * * * * * * * * * //

template<class Type>
void geometricFieldDST<Type>::operator+=
(
    const geometricFieldDST<Type>& gf
)
{
    GeometricField<Type, Foam::fvPatchField, Foam::volMesh>::operator+=(gf);
}

template<class Type>
void geometricFieldDST<Type>::operator+=
(
    const GeometricField<Type, Foam::fvPatchField, Foam::volMesh>& gf
)
{
    GeometricField<Type, Foam::fvPatchField, Foam::volMesh>::operator+=(gf);
}

template<class Type>
geometricFieldDST<Type>& geometricFieldDST<Type>::operator++(int)
{
    *this += (*deltaDST_);
    return *this;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace rhoDST






