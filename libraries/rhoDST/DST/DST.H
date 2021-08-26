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
                                      

\*---------------------------------------------------------------------------------------*/

#ifndef DST_H
#define DST_H

#include "fvMatrices.H"
#include "BlockLduSystem.H"


// These macros are taken from foam-extend-4.1
#define TEMPLATE template<template<class> class PatchField, class GeoMesh>

#define UNARY_FUNCTIONDST(Type1, Func)                                        \
                                                                              \
template<class Type1>                                                         \
void Func                                                                     \
(                                                                             \
    Field<Type1>& a,                                                          \
    const Field<Type1>& b                                                     \
);                                                                            \
                                                                              \
template<class Type1, template<class> class PatchField>                       \
void Func                                                                     \
(                                                                             \
    FieldField<PatchField, Type1>& a,                                         \
    const FieldField<PatchField, Type1>& b                                    \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<GeometricField<Type1, PatchField, GeoMesh> > Func                         \
(                                                                             \
    const GeometricField<Type1, PatchField, GeoMesh>& df                      \
);                                                                            \
                                                                              \
TEMPLATE                                                                      \
tmp<GeometricField<Type1, PatchField, GeoMesh> > Func                         \
(                                                                             \
    const tmp<GeometricField<Type1, PatchField, GeoMesh> >& tdf1              \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fvm functions Declaration
\*---------------------------------------------------------------------------*/

namespace DST
{    
    template<class vectorNType, class fieldType>
    void insert
    (
        Field<vectorNType>& vVf,
        const Field<fieldType>& vf,
        const direction& dir
    );
    
    template<class vectorNType, class fieldType>
    void insert
    (
        GeometricField<vectorNType, fvPatchField, volMesh>& vVf,
        const tmp<GeometricField<fieldType, fvPatchField, volMesh>>& vf,
        const direction& dir
    );

    template<class vectorNType, class fieldType>
    void insert
    (
        Field<vectorNType>& vVf,
        const tmp<Field<fieldType>>& vf,
        const direction& dir
    );
    
    template<class vectorNType, class fieldType>
    void insert
    (
        GeometricField<vectorNType, fvPatchField, volMesh>& vVf,
        const Field<fieldType>& vf,
        const direction& dir
    );

    template<class vectorNType, class fieldType>
    void retrieve
    (
        const GeometricField<vectorNType, fvPatchField, volMesh>& vVf,
        GeometricField<fieldType, fvPatchField, volMesh>& vf,
        const direction& dir
    );
    
    template<class jType, class fType>
    void insertScalarEqn
    (
        Field<jType>& J,
        const tmp<Field<fType>> tJi,
        int row,
        int col    
    );
    
    template<class jType, class fType>
    void insertScalarEqn
    (
        Field<jType>& J,
        const Field<fType>& Ji,
        int row,
        int col    
    );
    
    template<class jType, class fType>
    void insertVectorEqn
    (
        Field<jType>& J,
        const tmp<Field<fType>> tJi,
        int row,
        int col    
    );
    
    template<class jType, class fType>
    void insertVectorEqn
    (
        Field<jType>& J,
        const Field<fType>& Ji,
        int row,
        int col    
    );
    
    
    
    template<class Type>
    typename outerProduct< Type, Type >::type I();
    
    template<class Type>
    Type tensorI();
    
    inline scalar pos(const scalar s);

    inline scalar pos0(const scalar s);

    inline scalar neg(const scalar s);

    inline scalar neg0(const scalar s);
        
    UNARY_FUNCTIONDST(scalar, pos)
    UNARY_FUNCTIONDST(scalar, pos0)
    UNARY_FUNCTIONDST(scalar, neg)
    UNARY_FUNCTIONDST(scalar, neg0)
    
} // End namespace DST


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#undef TEMPLATE
#undef UNARY_FUNCTIONDST

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DST.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
