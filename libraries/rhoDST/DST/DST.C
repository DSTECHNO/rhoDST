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
                                      

\*---------------------------------------------------------------------------*/

#include "fvMesh.H"

// These macros are taken from foam-extend-4.1
#define TEMPLATE template<template<class> class PatchField, class GeoMesh>

#define UNARY_FUNCTIONDST(Type1, Func)                                        \
                                                                              \
template<class Type1>                                                         \
void Func                                                                     \
(                                                                             \
    Field<Type1>& a,                                                          \
    const Field<Type1>& b                                                     \
)                                                                             \
{                                                                             \
                                                                              \
    Type1* const __restrict__ fa = (a).begin();                               \
    const Type1* const __restrict__ fb = (b).begin();                         \
                                                                              \
        const label _n##i = (a).size();                                       \
        for (label i = 0; i < _n##i; i++)                                     \
        {                                                                     \
            fa[i] = DST::Func(fb[i]);                                         \
        }                                                                     \
}                                                                             \
                                                                              \
template<class Type1, template<class> class PatchField>                       \
void Func                                                                     \
(                                                                             \
    FieldField<PatchField, Type1>& a,                                         \
    const FieldField<PatchField, Type1>& b                                    \
)                                                                             \
{                                                                             \
    forAll(a,patchI)                                                          \
        DST::Func(a[patchI],b[patchI]);                                       \
}                                                                             \
                                                                              \
TEMPLATE                                                                      \
void Func                                                                     \
(                                                                             \
    GeometricField<Type1, PatchField, GeoMesh>& res,                          \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1                     \
)                                                                             \
{                                                                             \
    DST::Func(res.internalField(),gf1.internalField());                       \
    FieldField<PatchField, Type1>& a = res.boundaryField();                \
    const FieldField<PatchField, Type1>& b = gf1.boundaryField();          \
    DST::Func(a,b);                                                           \
}                                                                             \
                                                                              \
TEMPLATE                                                                      \
tmp<GeometricField<Type1, PatchField, GeoMesh> > Func                         \
(                                                                             \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1                     \
)                                                                             \
{                                                                             \
    tmp<GeometricField<Type1, PatchField, GeoMesh> > tRes                     \
    (                                                                         \
        new GeometricField<Type1, PatchField, GeoMesh>                        \
        (                                                                     \
            IOobject                                                          \
            (                                                                 \
                #Func "(" + gf1.name() + ')',                                 \
                gf1.instance(),                                               \
                gf1.db(),                                                     \
                IOobject::NO_READ,                                            \
                IOobject::NO_WRITE                                            \
            ),                                                                \
            gf1.mesh(),                                                       \
            dimless                                                           \
        )                                                                     \
    );                                                                        \
                                                                              \
    DST::Func(tRes(), gf1);                                                   \
                                                                              \
    return tRes;                                                              \
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace DST
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class vectorNType, class fieldType>
void insert
(
    Field<vectorNType>& vVfIn,
    const Field<fieldType>& vf,
    const direction& dir
)
{
    const direction nCmpts = pTraits<fieldType>::nComponents;
    direction localDir = dir;


    for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
    {
        scalarField vfInCurr(vf.component(cmptI));

        forAll (vfInCurr, cellI)
        {
            vVfIn[cellI](localDir) = vfInCurr[cellI];
        }

        localDir++;
    }
}

template<class vectorNType, class fieldType>
void insert
(
    GeometricField<vectorNType, fvPatchField, volMesh>& vVf,
    const tmp<GeometricField<fieldType, fvPatchField, volMesh>>& vf,
    const direction& dir
)
{
    insert(vVf.internalField(), vf().internalField(), dir);
    vf.clear();
}

template<class vectorNType, class fieldType>
void insert
(
    Field<vectorNType>& vVf,
    const tmp<Field<fieldType>>& vf,
    const direction& dir
)
{
    insert(vVf, vf(), dir);
    vf.clear();
}

template<class vectorNType, class fieldType>
void insert
(
    GeometricField<vectorNType, fvPatchField, volMesh>& vVf,
    const Field<fieldType>& vf,
    const direction& dir
)
{
    insert(vVf.internalField(), vf, dir);
}

template<class vectorNType, class fieldType>
void retrieve
(
    const GeometricField<vectorNType, fvPatchField, volMesh>& vVf,
    GeometricField<fieldType, fvPatchField, volMesh>& vf,
    const direction& dir
)
{
    const direction nCmpts = pTraits<fieldType>::nComponents;
    direction localDir = dir;

    const Field<vectorNType>& vVfIn = vVf.internalField();
    Field<fieldType>& vfIn = vf.internalField();

    for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
    {
        scalarField vfInCurr(vfIn.component(cmptI));

        forAll (vfInCurr, cellI)
        {
            vfInCurr[cellI] = vVfIn[cellI](localDir);
        }

        vfIn.replace(cmptI, vfInCurr);

        localDir++;
    }
}

template<class jType, class fType>
void insertScalarEqn
(
    Field<jType>& J,
    const tmp<Field<fType>> tJi,
    int row,
    int col    
)
{
    insertScalarEqn(J, tJi(), row, col);
}

template<class jType, class fType>
void insertScalarEqn
(
    Field<jType>& J,
    const Field<fType>& Ji,
    int row,
    int col    
)
{
    int nCols = pow(pTraits<jType>::nComponents,0.5);
    int nComponents = pTraits<fType>::nComponents;
    int dir = nCols*row + col;
    for(int i = 0; i<nComponents; i++)
    {
        J.replace(dir, Ji.component(i));
        dir++;
    }
}

template<class jType, class fType>
void insertVectorEqn
(
    Field<jType>& J,
    const tmp<Field<fType>> tJi,
    int row,
    int col    
)
{
    insertVectorEqn(J, tJi(), row, col);
}

template<class jType, class fType>
void insertVectorEqn
(
    Field<jType>& J,
    const Field<fType>& Ji,
    int row,
    int col    
)
{
    int nCols = pow(pTraits<jType>::nComponents,0.5);
    int dir = nCols*row + col;
    
    int rank = pTraits<fType>::rank;
    int rowF = 1;
    int colF = 1;
    int nComponents = pTraits<fType>::nComponents;
    int n = pow(nComponents,1./rank);
    
    if(rank == 1)
    {
        rowF = nComponents;
        colF = 1;
    }
    else if (rank == 2)
    {
        rowF = n;
        colF = n;
    }
    int nDif = nCols - colF;
    
    for(int i = 0; i<colF; i++)
    {
        for(int j =i; j <= colF*(rowF-1) + i; j += colF)
        {
            J.replace(dir, Ji.component(j));
            dir += nCols;
        }
        dir -= nCols*rowF -1;
    }
}

template<class Type>
typename outerProduct< Type, Type >::type I()
{
    typename outerProduct< Type, Type >::type a(outerProduct< Type, Type >::type::zero);
    int nComponents = pTraits<Type>::nComponents;
    for(int i=0; i<nComponents; i++)
        a(i,i) = 1;
    return a;
}

template<class Type>
Type tensorI()
{
    Type a(Type::zero);
    int nComponents = pow(pTraits<Type>::nComponents,0.5);
    for(int i=0; i<nComponents; i++)
        a(i,i) = 1;
    return a;
}

inline scalar pos(const scalar s)
{
    return (s > 0)? 1: 0;
}

inline scalar pos0(const scalar s)
{
    return (s >= 0)? 1: 0;
}

inline scalar neg(const scalar s)
{
    return (s < 0)? 1: 0;
}

inline scalar neg0(const scalar s)
{
    return (s <= 0)? 1: 0;
}

UNARY_FUNCTIONDST(scalar, pos)
UNARY_FUNCTIONDST(scalar, pos0)
UNARY_FUNCTIONDST(scalar, neg)
UNARY_FUNCTIONDST(scalar, neg0)
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace DST

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#undef TEMPLATE
#undef UNARY_FUNCTIONDST

// ************************************************************************* //
