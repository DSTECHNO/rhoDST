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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "Vector5typedefs.H"
#include "addToRunTimeSelectionTable.H"

#define forVector5Type(m, args...)              \
    m(vector5, Vector5, args)

#define forTensor5Type(m, args...)              \
    m(tensor5, Tensor5, args)

#define forDiagTensor5Type(m, args...)          \
    m(diagTensor5, DiagTensor5, args)

#define forSphericalTensor5Type(m, args...)     \
    m(sphericalTensor5, SphericalTensor5, args)

#define forVectorTensor5Types(m, args...)                            \
    m(tensor5, diagTensor5, sphericalTensor5, vector5, scalar, args)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

#define makeListType(type, Type, args...)                               \
    defineCompoundTypeName(List<type>, Type##List);                     \
    addCompoundToRunTimeSelectionTable(List<type>, type##List);

forVector5Type(makeListType)
forTensor5Type(makeListType)
forDiagTensor5Type(makeListType)
forSphericalTensor5Type(makeListType)

#undef makeListType

}

#define TEMPLATE
#include "FieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define VectorN_FieldFunctions(tensorType,diagTensorType,                    \
                               sphericalTensorType,vectorType,CmptType,      \
                               args...)                                      \
                                                                             \
UNARY_FUNCTION(CmptType, vectorType, cmptSum)                                \
                                                                             \
BINARY_FUNCTION(vectorType, vectorType, vectorType, cmptMultiply)            \
BINARY_TYPE_FUNCTION(vectorType, vectorType, vectorType, cmptMultiply)       \
                                                                             \
BINARY_OPERATOR(vectorType, CmptType, vectorType, /, divide)                 \
BINARY_TYPE_OPERATOR(vectorType, CmptType, vectorType, /, divide)            \
                                                                             \
BINARY_OPERATOR(vectorType, vectorType, vectorType, +, add)                  \
BINARY_OPERATOR(vectorType, vectorType, vectorType, -, subtract)             \
                                                                             \
BINARY_TYPE_OPERATOR(vectorType, vectorType, vectorType, +, add)             \
BINARY_TYPE_OPERATOR(vectorType, vectorType, vectorType, -, subtract)        \

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

forVectorTensor5Types(VectorN_FieldFunctions)

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef VectorN_FieldFunctions

#include "undefFieldFunctionsM.H"

#include "transformField.H"

#define TEMPLATE
#include "FieldFunctionsM.C"

#include "VectorN.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_TEMPLATE_FUNCTION(ReturnType, Type, Func)                      \
template<class Cmpt, int length>                                             \
inline void Func(ReturnType<Cmpt, length>&, const Type<Cmpt, length>&);

#define UNARY_TEMPLATE_FUNCTION_VS(ReturnType, Func)                         \
template<class Cmpt, int length>                                             \
inline void Func(ReturnType<Cmpt, length>&, const Cmpt&);

#define UNARY_TEMPLATE_FUNCTION_SV(Type, Func)                               \
template<class Cmpt, int length>                                             \
inline void Func(Cmpt&, const Type<Cmpt, length>&);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

UNARY_TEMPLATE_FUNCTION_SV(TensorN, contractScalar)
UNARY_TEMPLATE_FUNCTION_SV(DiagTensorN, contractScalar)
UNARY_TEMPLATE_FUNCTION_SV(SphericalTensorN, contractScalar)
UNARY_TEMPLATE_FUNCTION_SV(VectorN, contractScalar)

UNARY_TEMPLATE_FUNCTION(VectorN, TensorN, contractLinear)
UNARY_TEMPLATE_FUNCTION(VectorN, DiagTensorN, contractLinear)
UNARY_TEMPLATE_FUNCTION(VectorN, SphericalTensorN, contractLinear)

UNARY_TEMPLATE_FUNCTION_VS(VectorN, expandScalar)
UNARY_TEMPLATE_FUNCTION_VS(TensorN, expandScalar)
UNARY_TEMPLATE_FUNCTION_VS(DiagTensorN, expandScalar)
UNARY_TEMPLATE_FUNCTION_VS(SphericalTensorN, expandScalar)

UNARY_TEMPLATE_FUNCTION(TensorN, VectorN, expandLinear)
UNARY_TEMPLATE_FUNCTION(DiagTensorN, VectorN, expandLinear)
UNARY_TEMPLATE_FUNCTION(SphericalTensorN, VectorN, expandLinear)

UNARY_TEMPLATE_FUNCTION(VectorN, TensorN, sumToDiag)
UNARY_TEMPLATE_FUNCTION(VectorN, TensorN, sumMagToDiag)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ExpandTensorNI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef UNARY_TEMPLATE_FUNCTION
#undef UNARY_TEMPLATE_FUNCTION_VS
#undef UNARY_TEMPLATE_FUNCTION_SV

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define TensorN_FieldFunctions(tensorType,diagTensorType,                    \
                               sphericalTensorType,vectorType,CmptType,      \
                               args...)                                      \
                                                                             \
UNARY_FUNCTION(CmptType, tensorType, det)                                    \
UNARY_FUNCTION(tensorType, tensorType, inv)                                  \
UNARY_FUNCTION(diagTensorType, tensorType, diag)                             \
UNARY_FUNCTION(tensorType, tensorType, negSumDiag)                           \
                                                                             \
BINARY_OPERATOR(tensorType, CmptType, tensorType, /, divide)                 \
BINARY_TYPE_OPERATOR(tensorType, CmptType, tensorType, /, divide)            \
                                                                             \
BINARY_OPERATOR(vectorType, vectorType, tensorType, /, divide)               \
BINARY_TYPE_OPERATOR(vectorType, vectorType, tensorType, /, divide)          \
                                                                             \
BINARY_OPERATOR(tensorType, tensorType, tensorType, /, divide)               \
BINARY_TYPE_OPERATOR(tensorType, tensorType, tensorType, /, divide)          \
                                                                             \
BINARY_OPERATOR(tensorType, tensorType, diagTensorType, /, divide)           \
BINARY_TYPE_OPERATOR(tensorType, tensorType, diagTensorType, /, divide)      \
                                                                             \
BINARY_OPERATOR(tensorType, diagTensorType, tensorType, /, divide)           \
BINARY_TYPE_OPERATOR(tensorType, diagTensorType, tensorType, /, divide)      \
                                                                             \
BINARY_OPERATOR(tensorType, sphericalTensorType, tensorType, /, divide)      \
BINARY_TYPE_OPERATOR(tensorType, sphericalTensorType, tensorType, /, divide) \
                                                                             \
BINARY_OPERATOR(tensorType, tensorType, sphericalTensorType, /, divide)      \
BINARY_TYPE_OPERATOR(tensorType, tensorType, sphericalTensorType, /, divide) \
                                                                             \
BINARY_OPERATOR(tensorType, tensorType, tensorType, +, add)                  \
BINARY_OPERATOR(tensorType, tensorType, tensorType, -, subtract)             \
                                                                             \
BINARY_TYPE_OPERATOR(tensorType, tensorType, tensorType, +, add)             \
BINARY_TYPE_OPERATOR(tensorType, tensorType, tensorType, -, subtract)        \
                                                                             \
BINARY_OPERATOR(tensorType, diagTensorType, tensorType, +, add)              \
BINARY_OPERATOR(tensorType, diagTensorType, tensorType, -, subtract)         \
                                                                             \
BINARY_TYPE_OPERATOR(tensorType, diagTensorType, tensorType, +, add)         \
BINARY_TYPE_OPERATOR(tensorType, diagTensorType, tensorType, -, subtract)    \
                                                                             \
BINARY_OPERATOR(tensorType, tensorType, diagTensorType, +, add)              \
BINARY_OPERATOR(tensorType, tensorType, diagTensorType, -, subtract)         \
                                                                             \
BINARY_TYPE_OPERATOR(tensorType, tensorType, diagTensorType, +, add)         \
BINARY_TYPE_OPERATOR(tensorType, tensorType, diagTensorType, -, subtract)    \
                                                                             \
BINARY_OPERATOR(tensorType, sphericalTensorType, tensorType, +, add)         \
BINARY_OPERATOR(tensorType, sphericalTensorType, tensorType, -, subtract)    \
                                                                             \
BINARY_TYPE_OPERATOR(tensorType, sphericalTensorType, tensorType, +, add)    \
BINARY_TYPE_OPERATOR(tensorType, sphericalTensorType, tensorType, -, subtract) \
                                                                             \
BINARY_OPERATOR(tensorType, tensorType, sphericalTensorType, +, add)         \
BINARY_OPERATOR(tensorType, tensorType, sphericalTensorType, -, subtract)    \
                                                                             \
BINARY_TYPE_OPERATOR(tensorType, tensorType, sphericalTensorType, +, add)    \
BINARY_TYPE_OPERATOR(tensorType, tensorType, sphericalTensorType, -, subtract) \
                                                                             \
template<>                                                                   \
tmp<Field<tensorType> > transformFieldMask<tensorType>                       \
(                                                                            \
    const Field<diagTensorType>& dtf                                         \
)                                                                            \
{                                                                            \
    tmp<Field<tensorType> > tRes(new Field<tensorType>(dtf.size()));         \
    Field<tensorType>& res = tRes();                                         \
    TFOR_ALL_F_OP_F(tensorType, res, =, diagTensorType, dtf)                 \
    return tRes;                                                             \
}                                                                            \
                                                                             \
template<>                                                                   \
tmp<Field<tensorType> > transformFieldMask<tensorType>                       \
(                                                                            \
    const Field<sphericalTensorType>& stf                                    \
)                                                                            \
{                                                                            \
    tmp<Field<tensorType> > tRes(new Field<tensorType>(stf.size()));         \
    Field<tensorType>& res = tRes();                                         \
    TFOR_ALL_F_OP_F(tensorType, res, =, sphericalTensorType, stf)            \
    return tRes;                                                             \
}                                                                            \


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

forVectorTensor5Types(TensorN_FieldFunctions)

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef TensorN_FieldFunctions

#include "undefFieldFunctionsM.H"

#include "Field.H"
#include "FieldM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_FUNCTION(typeF1, typeF2, FUNC)                                 \
                                                                             \
void FUNC(Field<typeF1>& f1, const UList<typeF2>& f2)                        \
{                                                                            \
    checkFields(f1, f2, #FUNC "(f1,f2)");                                    \
                                                                             \
    List_ACCESS(typeF1, f1, f1P);                                            \
    List_CONST_ACCESS(typeF2, f2, f2P);                                      \
                                                                             \
    List_FOR_ALL(f1,i)                                                       \
        FUNC(List_ELEM(f1, f1P, i), List_ELEM(f2, f2P, i));                  \
    List_END_FOR_ALL                                                         \
}                                                                            \
                                                                             \
void FUNC(Field<typeF1>& f1, const tmp<Field<typeF2> >& tf2)                 \
{                                                                            \
     FUNC(f1,tf2());                                                         \
     tf2.clear();                                                            \
}

#define ExpandFieldFunctions(tensorType, diagTensorType, sphericalTensorType, \
        vectorType, cmptType, args...)                                        \
                                                                              \
UNARY_FUNCTION(cmptType, tensorType, contractScalar)                          \
UNARY_FUNCTION(cmptType, diagTensorType, contractScalar)                      \
UNARY_FUNCTION(cmptType, sphericalTensorType, contractScalar)                 \
UNARY_FUNCTION(cmptType, vectorType, contractScalar)                          \
                                                                              \
UNARY_FUNCTION(vectorType, tensorType, contractLinear)                        \
UNARY_FUNCTION(vectorType, diagTensorType, contractLinear)                    \
UNARY_FUNCTION(vectorType, sphericalTensorType, contractLinear)               \
                                                                              \
UNARY_FUNCTION(vectorType, cmptType, expandScalar)                            \
UNARY_FUNCTION(tensorType, cmptType, expandScalar)                            \
UNARY_FUNCTION(diagTensorType, cmptType, expandScalar)                        \
UNARY_FUNCTION(sphericalTensorType, cmptType, expandScalar)                   \
                                                                              \
UNARY_FUNCTION(tensorType, vectorType, expandLinear)                          \
UNARY_FUNCTION(diagTensorType, vectorType, expandLinear)                      \
UNARY_FUNCTION(sphericalTensorType, vectorType, expandLinear)                 \
                                                                              \
UNARY_FUNCTION(vectorType, tensorType, sumToDiag)                             \
UNARY_FUNCTION(vectorType, tensorType, sumMagToDiag)


namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

forVectorTensor5Types(ExpandFieldFunctions)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef UNARY_FUNCTION
#undef ExpandFieldFunctions

#define TEMPLATE
#include "FieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define DiagTensorN_FieldFunctions(tensorType, diagTensorType, sphericalTensorType, \
    vectorType, CmptType, args...)                                                  \
                                                                                    \
UNARY_FUNCTION(diagTensorType, diagTensorType, inv)                                 \
UNARY_FUNCTION(diagTensorType, diagTensorType, diag)                                \
                                                                                    \
BINARY_OPERATOR(diagTensorType, CmptType, diagTensorType, /, divide)                \
BINARY_TYPE_OPERATOR(diagTensorType, CmptType, diagTensorType, /, divide)           \
                                                                                    \
BINARY_OPERATOR(vectorType, vectorType, diagTensorType, /, divide)                  \
BINARY_TYPE_OPERATOR(vectorType, vectorType, diagTensorType, /, divide)             \
                                                                                    \
BINARY_OPERATOR(diagTensorType, diagTensorType, diagTensorType, /, divide)          \
BINARY_TYPE_OPERATOR(diagTensorType, diagTensorType, diagTensorType, /, divide)     \
                                                                                    \
BINARY_OPERATOR(diagTensorType, sphericalTensorType, diagTensorType, /, divide)     \
BINARY_TYPE_OPERATOR(diagTensorType, sphericalTensorType, diagTensorType, /, divide)\
                                                                                    \
BINARY_OPERATOR(diagTensorType, diagTensorType, sphericalTensorType, /, divide)     \
BINARY_TYPE_OPERATOR(diagTensorType, diagTensorType, sphericalTensorType, /, divide)\
                                                                                    \
BINARY_OPERATOR(diagTensorType, diagTensorType, diagTensorType, +, add)             \
BINARY_OPERATOR(diagTensorType, diagTensorType, diagTensorType, -, subtract)        \
                                                                                    \
BINARY_TYPE_OPERATOR(diagTensorType, diagTensorType, diagTensorType, +, add)        \
BINARY_TYPE_OPERATOR(diagTensorType, diagTensorType, diagTensorType, -, subtract)   \
                                                                                    \
BINARY_OPERATOR(diagTensorType, sphericalTensorType, diagTensorType, +, add)        \
BINARY_OPERATOR(diagTensorType, sphericalTensorType, diagTensorType, -, subtract)   \
                                                                                    \
BINARY_TYPE_OPERATOR(diagTensorType, sphericalTensorType, diagTensorType, +, add)   \
BINARY_TYPE_OPERATOR(diagTensorType, sphericalTensorType, diagTensorType, -, subtract)  \
                                                                                    \
BINARY_OPERATOR(diagTensorType, diagTensorType, sphericalTensorType, +, add)        \
BINARY_OPERATOR(diagTensorType, diagTensorType, sphericalTensorType, -, subtract)   \
                                                                                    \
BINARY_TYPE_OPERATOR(diagTensorType, diagTensorType, sphericalTensorType, +, add)   \
BINARY_TYPE_OPERATOR(diagTensorType, diagTensorType, sphericalTensorType, -, subtract)  \
                                                                                    \
template<>                                                                          \
tmp<Field<diagTensorType> > transformFieldMask<diagTensorType>                      \
(                                                                                   \
    const Field<sphericalTensorType>& stf                                           \
)                                                                                   \
{                                                                                   \
    tmp<Field<diagTensorType> > tRes( new Field<diagTensorType>(stf.size()) );      \
    Field<diagTensorType>& res = tRes();                                            \
    TFOR_ALL_F_OP_F(diagTensorType, res, =, sphericalTensorType, stf)               \
    return tRes;                                                                    \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

forVectorTensor5Types(DiagTensorN_FieldFunctions)

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef DiagTensorN_FieldFunctions

#include "undefFieldFunctionsM.H"


#define TEMPLATE
#include "FieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define SphericalTensorN_FieldFunctions(tensorType, diagTensorType,                 \
    sphericalTensorType, vectorType, CmptType, args...)                             \
                                                                                    \
UNARY_FUNCTION(sphericalTensorType, sphericalTensorType, inv)                       \
UNARY_FUNCTION(sphericalTensorType, sphericalTensorType, diag)                      \
                                                                                    \
BINARY_OPERATOR(sphericalTensorType, CmptType, sphericalTensorType, /, divide)      \
BINARY_TYPE_OPERATOR(sphericalTensorType, CmptType, sphericalTensorType, /, divide) \
                                                                                    \
BINARY_OPERATOR(vectorType, vectorType, sphericalTensorType, /, divide)             \
BINARY_TYPE_OPERATOR(vectorType, vectorType, sphericalTensorType, /, divide)        \
                                                                                    \
BINARY_OPERATOR(sphericalTensorType, sphericalTensorType, sphericalTensorType, /, divide)           \
BINARY_TYPE_OPERATOR(sphericalTensorType, sphericalTensorType, sphericalTensorType, /, divide)      \
                                                                                    \
BINARY_OPERATOR(sphericalTensorType, sphericalTensorType, sphericalTensorType, +, add)             \
BINARY_OPERATOR(sphericalTensorType, sphericalTensorType, sphericalTensorType, -, subtract)        \
                                                                                    \
BINARY_TYPE_OPERATOR(sphericalTensorType, sphericalTensorType, sphericalTensorType, +, add)        \
BINARY_TYPE_OPERATOR(sphericalTensorType, sphericalTensorType, sphericalTensorType, -, subtract)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

forVectorTensor5Types(SphericalTensorN_FieldFunctions)

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef SphericalTensorN_FieldFunctions

#include "undefFieldFunctionsM.H"

#define TEMPLATE template<template<class> class Field>
#include "FieldFieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define VectorN_FieldFunctions(tensorType, diagTensorType, sphericalTensorType,         \
    vectorType, CmptType, args...)                                                      \
                                                                                        \
UNARY_FUNCTION(CmptType, vectorType, cmptSum)                                           \
                                                                                        \
BINARY_FUNCTION(vectorType, vectorType, vectorType, cmptMultiply)                       \
BINARY_TYPE_FUNCTION(vectorType, vectorType, vectorType, cmptMultiply)                  \
                                                                                        \
BINARY_OPERATOR(vectorType, vectorType, CmptType, /,divide)                             \
BINARY_TYPE_OPERATOR(vectorType, vectorType, CmptType, /,divide)                        \
                                                                                        \
BINARY_OPERATOR(vectorType, vectorType, vectorType, +, add)                             \
BINARY_OPERATOR(vectorType, vectorType, vectorType, -, subtract)                        \
                                                                                        \
BINARY_TYPE_OPERATOR(vectorType, vectorType, vectorType, +, add)                        \
BINARY_TYPE_OPERATOR(vectorType, vectorType, vectorType, -, subtract)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

forVectorTensor5Types(VectorN_FieldFunctions)

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef VectorN_FieldFunctions

#include "undefFieldFunctionsM.H"

#define TEMPLATE template<template<class> class Field>
#include "FieldFieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define TensorN_FieldFunctions(tensorType, diagTensorType, sphericalTensorType,         \
    vectorType, CmptType, args...)                                                      \
                                                                                        \
UNARY_FUNCTION(tensorType, tensorType, inv)                                             \
UNARY_FUNCTION(diagTensorType, tensorType, diag)                                        \
UNARY_FUNCTION(tensorType, tensorType, negSumDiag)                                      \
UNARY_FUNCTION(vectorType, tensorType, contractLinear)                                  \
UNARY_FUNCTION(CmptType, tensorType, contractScalar)                                    \
                                                                                        \
BINARY_OPERATOR(tensorType, CmptType, tensorType, /, divide)                            \
BINARY_TYPE_OPERATOR(tensorType, CmptType, tensorType, /, divide)                       \
                                                                                        \
BINARY_OPERATOR(vectorType, vectorType, tensorType, /, divide)                          \
BINARY_TYPE_OPERATOR(vectorType, vectorType, tensorType, /, divide)                     \
                                                                                        \
BINARY_OPERATOR(tensorType, tensorType, tensorType, /, divide)                          \
BINARY_TYPE_OPERATOR(tensorType, tensorType, tensorType, /, divide)                     \
                                                                                        \
BINARY_OPERATOR(tensorType, tensorType, diagTensorType, /, divide)                      \
BINARY_TYPE_OPERATOR(tensorType, tensorType, diagTensorType, /, divide)                 \
                                                                                        \
BINARY_OPERATOR(tensorType, diagTensorType, tensorType, /, divide)                      \
BINARY_TYPE_OPERATOR(tensorType, diagTensorType, tensorType, /, divide)                 \
                                                                                        \
BINARY_OPERATOR(tensorType, sphericalTensorType, tensorType, /, divide)                 \
BINARY_TYPE_OPERATOR(tensorType, sphericalTensorType, tensorType, /, divide)            \
                                                                                        \
BINARY_OPERATOR(tensorType, tensorType, sphericalTensorType, /, divide)                 \
BINARY_TYPE_OPERATOR(tensorType, tensorType, sphericalTensorType, /, divide)            \
                                                                                        \
BINARY_OPERATOR(tensorType, tensorType, tensorType, +, add)                             \
BINARY_OPERATOR(tensorType, tensorType, tensorType, -, subtract)                        \
                                                                                        \
BINARY_TYPE_OPERATOR(tensorType, tensorType, tensorType, +, add)                        \
BINARY_TYPE_OPERATOR(tensorType, tensorType, tensorType, -, subtract)                   \
                                                                                        \
BINARY_OPERATOR(tensorType, diagTensorType, tensorType, +, add)                         \
BINARY_OPERATOR(tensorType, diagTensorType, tensorType, -, subtract)                    \
                                                                                        \
BINARY_TYPE_OPERATOR(tensorType, diagTensorType, tensorType, +, add)                    \
BINARY_TYPE_OPERATOR(tensorType, diagTensorType, tensorType, -, subtract)               \
                                                                                        \
BINARY_OPERATOR(tensorType, sphericalTensorType, tensorType, +, add)                    \
BINARY_OPERATOR(tensorType, sphericalTensorType, tensorType, -, subtract)               \
                                                                                        \
BINARY_TYPE_OPERATOR(tensorType, sphericalTensorType, tensorType, +, add)               \
BINARY_TYPE_OPERATOR(tensorType, sphericalTensorType, tensorType, -, subtract)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

forVectorTensor5Types(TensorN_FieldFunctions)

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef TensorN_FieldFunctions

#include "undefFieldFunctionsM.H"

#define TEMPLATE template<template<class> class Field>
#include "FieldFieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define DiagTensorN_FieldFunctions(tensorType, diagTensorType, sphericalTensorType, \
    vectorType, CmptType, args...)                                                  \
                                                                                    \
UNARY_FUNCTION(diagTensorType, diagTensorType, inv)                                 \
UNARY_FUNCTION(diagTensorType, diagTensorType, diag)                                \
UNARY_FUNCTION(vectorType, diagTensorType, contractLinear)                          \
UNARY_FUNCTION(CmptType, diagTensorType, contractScalar)                            \
                                                                                    \
BINARY_OPERATOR(diagTensorType, CmptType, diagTensorType, /, divide)                \
BINARY_TYPE_OPERATOR(diagTensorType, CmptType, diagTensorType, /, divide)           \
                                                                                    \
BINARY_OPERATOR(vectorType, vectorType, diagTensorType, /, divide)                  \
BINARY_TYPE_OPERATOR(vectorType, vectorType, diagTensorType, /, divide)             \
                                                                                    \
BINARY_OPERATOR(diagTensorType, diagTensorType, diagTensorType, /, divide)          \
BINARY_TYPE_OPERATOR(diagTensorType, diagTensorType, diagTensorType, /, divide)     \
                                                                                    \
BINARY_OPERATOR(diagTensorType, sphericalTensorType, diagTensorType, /, divide)     \
BINARY_TYPE_OPERATOR(diagTensorType, sphericalTensorType, diagTensorType, /, divide)\
                                                                                    \
BINARY_OPERATOR(diagTensorType, diagTensorType, sphericalTensorType, /, divide)     \
BINARY_TYPE_OPERATOR(diagTensorType, diagTensorType, sphericalTensorType, /, divide)\
                                                                                    \
BINARY_OPERATOR(diagTensorType, diagTensorType, diagTensorType, +, add)             \
BINARY_OPERATOR(diagTensorType, diagTensorType, diagTensorType, -, subtract)        \
                                                                                    \
BINARY_TYPE_OPERATOR(diagTensorType, diagTensorType, diagTensorType, +, add)        \
BINARY_TYPE_OPERATOR(diagTensorType, diagTensorType, diagTensorType, -, subtract)   \
                                                                                    \
BINARY_OPERATOR(diagTensorType, sphericalTensorType, diagTensorType, +, add)        \
BINARY_OPERATOR(diagTensorType, sphericalTensorType, diagTensorType, -, subtract)   \
                                                                                    \
BINARY_TYPE_OPERATOR(diagTensorType, sphericalTensorType, diagTensorType, +, add)   \
BINARY_TYPE_OPERATOR(diagTensorType, sphericalTensorType, diagTensorType, -, subtract)  \
                                                                                    \
BINARY_OPERATOR(diagTensorType, diagTensorType, sphericalTensorType, +, add)        \
BINARY_OPERATOR(diagTensorType, diagTensorType, sphericalTensorType, -, subtract)   \
                                                                                    \
BINARY_TYPE_OPERATOR(diagTensorType, diagTensorType, sphericalTensorType, +, add)   \
BINARY_TYPE_OPERATOR(diagTensorType, diagTensorType, sphericalTensorType, -, subtract)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

forVectorTensor5Types(DiagTensorN_FieldFunctions)

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef DiagTensorN_FieldFunctions

#include "undefFieldFunctionsM.H"

#define TEMPLATE template<template<class> class Field>
#include "FieldFieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define SphericalTensorN_FieldFunctions(tensorType, diagTensorType,                 \
    sphericalTensorType, vectorType, CmptType, args...)                             \
                                                                                    \
UNARY_FUNCTION(sphericalTensorType, sphericalTensorType, inv)                       \
UNARY_FUNCTION(sphericalTensorType, sphericalTensorType, diag)                      \
UNARY_FUNCTION(CmptType, sphericalTensorType, contractLinear)                       \
UNARY_FUNCTION(CmptType, sphericalTensorType, contractScalar)                       \
                                                                                    \
BINARY_OPERATOR(sphericalTensorType, CmptType, sphericalTensorType, /, divide)      \
BINARY_TYPE_OPERATOR(sphericalTensorType, CmptType, sphericalTensorType, /, divide) \
                                                                                    \
BINARY_OPERATOR(vectorType, vectorType, sphericalTensorType, /, divide)             \
BINARY_TYPE_OPERATOR(vectorType, vectorType, sphericalTensorType, /, divide)        \
                                                                                    \
BINARY_OPERATOR(sphericalTensorType, sphericalTensorType, sphericalTensorType, /, divide)           \
BINARY_TYPE_OPERATOR(sphericalTensorType, sphericalTensorType, sphericalTensorType, /, divide)      \
                                                                                    \
BINARY_OPERATOR(sphericalTensorType, sphericalTensorType, sphericalTensorType, +, add)             \
BINARY_OPERATOR(sphericalTensorType, sphericalTensorType, sphericalTensorType, -, subtract)        \
                                                                                    \
BINARY_TYPE_OPERATOR(sphericalTensorType, sphericalTensorType, sphericalTensorType, +, add)        \
BINARY_TYPE_OPERATOR(sphericalTensorType, sphericalTensorType, sphericalTensorType, -, subtract)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

forVectorTensor5Types(SphericalTensorN_FieldFunctions)

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef SphericalTensorN_FieldFunctions

#include "undefFieldFunctionsM.H"

#define TEMPLATE template<class GeoMesh>
#include "DimensionedFieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define VectorN_FieldFunctions(tensorType, diagTensorType, sphericalTensorType, \
    vectorType, CmptType, args...)                                              \
                                                                                \
UNARY_FUNCTION(CmptType, vectorType, cmptSum, cmptSum)                          \
                                                                                \
BINARY_FUNCTION(vectorType, vectorType, vectorType, cmptMultiply)               \
BINARY_TYPE_FUNCTION(vectorType, vectorType, vectorType, cmptMultiply)          \
                                                                                \
BINARY_OPERATOR(vectorType, CmptType, vectorType, /,'|',divide)                 \
BINARY_TYPE_OPERATOR(vectorType, CmptType, vectorType, /,'|',divide)            \
                                                                                \
BINARY_OPERATOR(vectorType, vectorType, vectorType, +,'+',add)                  \
BINARY_OPERATOR(vectorType, vectorType, vectorType, -,'-',subtract)             \
                                                                                \
BINARY_TYPE_OPERATOR(vectorType, vectorType, vectorType, +,'+',add)             \
BINARY_TYPE_OPERATOR(vectorType, vectorType, vectorType, -,'-',subtract)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

forVectorTensor5Types(VectorN_FieldFunctions)

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef VectorN_FieldFunctions
#include "undefFieldFunctionsM.H"

#define TEMPLATE template<class GeoMesh>
#include "DimensionedFieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define TensorN_FieldFunctions(tensorType, diagTensorType, sphericalTensorType,     \
    vectorType, CmptType, args...)                                                  \
                                                                                    \
UNARY_FUNCTION(tensorType, tensorType,inv,inv)                                      \
UNARY_FUNCTION(diagTensorType, tensorType,diag,diag)                                \
UNARY_FUNCTION(tensorType, tensorType, negSumDiag, negSumDiag)                      \
UNARY_FUNCTION(vectorType, tensorType, contractLinear,contractLinear)               \
UNARY_FUNCTION(CmptType, tensorType, contractScalar,contractLinear)                 \
                                                                                    \
BINARY_OPERATOR(tensorType, CmptType, tensorType, /,'|',divide)                     \
BINARY_TYPE_OPERATOR(tensorType, CmptType, tensorType, /,'|',divide)                \
                                                                                    \
BINARY_OPERATOR(vectorType, vectorType, tensorType, /,'|',divide)                   \
BINARY_TYPE_OPERATOR(vectorType, vectorType, tensorType, /,'|',divide)              \
                                                                                    \
BINARY_OPERATOR(tensorType, tensorType, tensorType, /,'|',divide)                   \
BINARY_TYPE_OPERATOR(tensorType, tensorType, tensorType, /,'|',divide)              \
                                                                                    \
BINARY_OPERATOR(tensorType, tensorType, diagTensorType, /,'|',divide)               \
BINARY_TYPE_OPERATOR(tensorType, tensorType, diagTensorType, /,'|',divide)          \
                                                                                    \
BINARY_OPERATOR(tensorType, diagTensorType, tensorType, /,'|',divide)               \
BINARY_TYPE_OPERATOR(tensorType, diagTensorType, tensorType, /,'|',divide)          \
                                                                                    \
BINARY_OPERATOR(tensorType, sphericalTensorType, tensorType, /,'|',divide)          \
BINARY_TYPE_OPERATOR(tensorType, sphericalTensorType, tensorType, /,'|',divide)     \
                                                                                    \
BINARY_OPERATOR(tensorType, tensorType, sphericalTensorType, /,'|',divide)          \
BINARY_TYPE_OPERATOR(tensorType, tensorType, sphericalTensorType, /,'|',divide)     \
                                                                                    \
BINARY_OPERATOR(tensorType, tensorType, tensorType, +,'+',add)                      \
BINARY_OPERATOR(tensorType, tensorType, tensorType, -,'-',subtract)                 \
                                                                                    \
BINARY_TYPE_OPERATOR(tensorType, tensorType, tensorType, +,'+',add)                 \
BINARY_TYPE_OPERATOR(tensorType, tensorType, tensorType, -,'-',subtract)            \
                                                                                    \
BINARY_OPERATOR(tensorType, diagTensorType, tensorType, +,'+',add)                  \
BINARY_OPERATOR(tensorType, diagTensorType, tensorType, -,'-',subtract)             \
                                                                                    \
BINARY_TYPE_OPERATOR(tensorType, diagTensorType, tensorType, +,'+',add)             \
BINARY_TYPE_OPERATOR(tensorType, diagTensorType, tensorType, -,'-',subtract)        \
                                                                                    \
BINARY_OPERATOR(tensorType, sphericalTensorType, tensorType, +,'+',add)             \
BINARY_OPERATOR(tensorType, sphericalTensorType, tensorType, -,'-',subtract)        \
                                                                                    \
BINARY_TYPE_OPERATOR(tensorType, sphericalTensorType, tensorType, +,'+',add)        \
BINARY_TYPE_OPERATOR(tensorType, sphericalTensorType, tensorType, -,'-',subtract)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

forVectorTensor5Types(TensorN_FieldFunctions)

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef TensorN_FieldFunctions

#include "undefFieldFunctionsM.H"

#define TEMPLATE template<class GeoMesh>
#include "DimensionedFieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define DiagTensorN_FieldFunctions(tensorType, diagTensorType, sphericalTensorType,     \
    vectorType, CmptType, args...)                                                      \
                                                                                        \
UNARY_FUNCTION(diagTensorType, diagTensorType,inv,inv)                                  \
UNARY_FUNCTION(diagTensorType, diagTensorType,diag,diag)                                \
UNARY_FUNCTION(vectorType, diagTensorType, contractLinear,contractLinear)               \
UNARY_FUNCTION(CmptType, diagTensorType, contractScalar,contractLinear)                 \
                                                                                        \
BINARY_OPERATOR(diagTensorType, CmptType, diagTensorType, /,'|',divide)                 \
BINARY_TYPE_OPERATOR(diagTensorType, CmptType, diagTensorType, /,'|',divide)            \
                                                                                        \
BINARY_OPERATOR(vectorType, vectorType, diagTensorType, /,'|',divide)                   \
BINARY_TYPE_OPERATOR(vectorType, vectorType, diagTensorType, /,'|',divide)              \
                                                                                        \
BINARY_OPERATOR(diagTensorType, diagTensorType, diagTensorType, /,'|',divide)           \
BINARY_TYPE_OPERATOR(diagTensorType, diagTensorType, diagTensorType, /,'|',divide)      \
                                                                                        \
BINARY_OPERATOR(diagTensorType, sphericalTensorType, diagTensorType, /,'|',divide)      \
BINARY_TYPE_OPERATOR(diagTensorType, sphericalTensorType, diagTensorType, /,'|',divide) \
                                                                                        \
BINARY_OPERATOR(diagTensorType, diagTensorType, sphericalTensorType, /,'|',divide)      \
BINARY_TYPE_OPERATOR(diagTensorType, diagTensorType, sphericalTensorType, /,'|',divide) \
                                                                                        \
BINARY_OPERATOR(diagTensorType, diagTensorType, diagTensorType, +,'+',add)              \
BINARY_OPERATOR(diagTensorType, diagTensorType, diagTensorType, -,'-',subtract)         \
                                                                                        \
BINARY_TYPE_OPERATOR(diagTensorType, diagTensorType, diagTensorType, +,'+', add)        \
BINARY_TYPE_OPERATOR(diagTensorType, diagTensorType, diagTensorType, -,'-', subtract)   \
                                                                                        \
BINARY_OPERATOR(diagTensorType, sphericalTensorType, diagTensorType, +,'+', add)        \
BINARY_OPERATOR(diagTensorType, sphericalTensorType, diagTensorType, -,'-', subtract)   \
                                                                                        \
BINARY_TYPE_OPERATOR(diagTensorType, sphericalTensorType, diagTensorType, +,'+', add)   \
BINARY_TYPE_OPERATOR(diagTensorType, sphericalTensorType, diagTensorType, -,'-', subtract)  \
                                                                                        \
BINARY_OPERATOR(diagTensorType, diagTensorType, sphericalTensorType, +,'+', add)        \
BINARY_OPERATOR(diagTensorType, diagTensorType, sphericalTensorType, -,'-', subtract)   \
                                                                                        \
BINARY_TYPE_OPERATOR(diagTensorType, diagTensorType, sphericalTensorType, +,'+', add)   \
BINARY_TYPE_OPERATOR(diagTensorType, diagTensorType, sphericalTensorType, -,'+', subtract)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

forVectorTensor5Types(DiagTensorN_FieldFunctions)

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef DiagTensorN_FieldFunctions

#include "undefFieldFunctionsM.H"

#define TEMPLATE template<class GeoMesh>
#include "DimensionedFieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define SphericalTensorN_FieldFunctions(tensorType, diagTensorType,                 \
    sphericalTensorType, vectorType, CmptType, args...)                             \
                                                                                    \
UNARY_FUNCTION(sphericalTensorType, sphericalTensorType,inv,inv)                    \
UNARY_FUNCTION(sphericalTensorType, sphericalTensorType,diag,diag)                  \
UNARY_FUNCTION(CmptType, sphericalTensorType, contractLinear,contractLinear)        \
UNARY_FUNCTION(CmptType, sphericalTensorType, contractScalar,contractLinear)        \
                                                                                    \
BINARY_OPERATOR(sphericalTensorType, CmptType, sphericalTensorType, /,'|',divide)   \
BINARY_TYPE_OPERATOR(sphericalTensorType, CmptType, sphericalTensorType, /,'|',divide) \
                                                                                    \
BINARY_OPERATOR(vectorType, vectorType, sphericalTensorType, /,'|',divide)          \
BINARY_TYPE_OPERATOR(vectorType, vectorType, sphericalTensorType, /,'|',divide)     \
                                                                                    \
BINARY_OPERATOR(sphericalTensorType, sphericalTensorType, sphericalTensorType, /,'|',divide)           \
BINARY_TYPE_OPERATOR(sphericalTensorType, sphericalTensorType, sphericalTensorType, /,'|',divide)      \
                                                                                    \
BINARY_OPERATOR(sphericalTensorType, sphericalTensorType, sphericalTensorType, +,'+',add)             \
BINARY_OPERATOR(sphericalTensorType, sphericalTensorType, sphericalTensorType, -,'-',subtract)        \
                                                                                    \
BINARY_TYPE_OPERATOR(sphericalTensorType, sphericalTensorType, sphericalTensorType, +,'+',add)        \
BINARY_TYPE_OPERATOR(sphericalTensorType, sphericalTensorType, sphericalTensorType, -,'-',subtract)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

forVectorTensor5Types(SphericalTensorN_FieldFunctions)

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef SphericalTensorN_FieldFunctions

#include "undefFieldFunctionsM.H"

#define TEMPLATE template<template<class> class PatchField, class GeoMesh>
#include "GeometricFieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define VectorN_FieldFunctions(tensorType, diagTensorType, sphericalTensorType, \
    vectorType, CmptType, args...)                                              \
                                                                                \
UNARY_FUNCTION(CmptType, vectorType, cmptSum, cmptSum)                          \
                                                                                \
BINARY_FUNCTION(vectorType, vectorType, vectorType, cmptMultiply)               \
BINARY_TYPE_FUNCTION(vectorType, vectorType, vectorType, cmptMultiply)          \
                                                                                \
BINARY_OPERATOR(vectorType, CmptType, vectorType, *,'*',multiply)               \
BINARY_TYPE_OPERATOR(vectorType, CmptType, vectorType, *,'*',multiply)          \
                                                                                \
BINARY_OPERATOR(vectorType, CmptType, vectorType, /,'|',divide)                 \
BINARY_TYPE_OPERATOR(vectorType, CmptType, vectorType, /,'|',divide)            \
                                                                                \
BINARY_OPERATOR(vectorType, vectorType, vectorType, +,'+',add)                  \
BINARY_OPERATOR(vectorType, vectorType, vectorType, -,'-',subtract)             \
                                                                                \
BINARY_TYPE_OPERATOR(vectorType, vectorType, vectorType, +,'+',add)             \
BINARY_TYPE_OPERATOR(vectorType, vectorType, vectorType, -,'-',subtract)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

forVectorTensor5Types(VectorN_FieldFunctions)

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef VectorN_FieldFunctions

#include "undefFieldFunctionsM.H"

#include "Field.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_FUNCTION(typeF1, typeF2, FUNC)                                  \
void FUNC(Field<typeF1>& f1, const UList<typeF2>& f2);                        \
void FUNC(Field<typeF1>& f1, const tmp<Field<typeF2> >& tf2);


#define ExpandFieldFunctions(tensorType, diagTensorType, sphericalTensorType, \
        vectorType, cmptType, args...)                                        \
                                                                              \
UNARY_FUNCTION(cmptType, tensorType, contractScalar)                          \
UNARY_FUNCTION(cmptType, diagTensorType, contractScalar)                      \
UNARY_FUNCTION(cmptType, sphericalTensorType, contractScalar)                 \
UNARY_FUNCTION(cmptType, vectorType, contractScalar)                          \
                                                                              \
UNARY_FUNCTION(vectorType, tensorType, contractLinear)                        \
UNARY_FUNCTION(vectorType, diagTensorType, contractLinear)                    \
UNARY_FUNCTION(vectorType, sphericalTensorType, contractLinear)               \
                                                                              \
UNARY_FUNCTION(vectorType, cmptType, expandScalar)                            \
UNARY_FUNCTION(tensorType, cmptType, expandScalar)                            \
UNARY_FUNCTION(diagTensorType, cmptType, expandScalar)                        \
UNARY_FUNCTION(sphericalTensorType, cmptType, expandScalar)                   \
                                                                              \
UNARY_FUNCTION(tensorType, vectorType, expandLinear)                          \
UNARY_FUNCTION(diagTensorType, vectorType, expandLinear)                      \
UNARY_FUNCTION(sphericalTensorType, vectorType, expandLinear)                 \
                                                                              \
UNARY_FUNCTION(vectorType, tensorType, sumToDiag)                             \
UNARY_FUNCTION(vectorType, tensorType, sumMagToDiag)


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

namespace Foam
{

forVectorTensor5Types(ExpandFieldFunctions)

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef UNARY_FUNCTION
#undef ExpandFieldFunctions

#define TEMPLATE template<template<class> class PatchField, class GeoMesh>
#include "GeometricFieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define TensorN_FieldFunctions(tensorType, diagTensorType, sphericalTensorType,     \
    vectorType, CmptType, args...)                                                  \
                                                                                    \
UNARY_FUNCTION(tensorType, tensorType, T, transform)                                \
UNARY_FUNCTION(diagTensorType, tensorType,diag,diag)                                \
UNARY_FUNCTION(tensorType, tensorType,negSumDiag,negSumDiag)                        \
UNARY_FUNCTION(CmptType, tensorType,contractScalar,contractScalar)                  \
UNARY_FUNCTION(vectorType, tensorType,contractLinear,contractLinear)                \
                                                                                    \
BINARY_OPERATOR(tensorType, CmptType, tensorType, *,'*',multiply)                   \
BINARY_TYPE_OPERATOR(tensorType, CmptType, tensorType, *,'*',multiply)              \
                                                                                    \
BINARY_OPERATOR(tensorType, CmptType, tensorType, /,'|',divide)                     \
BINARY_TYPE_OPERATOR(tensorType, CmptType, tensorType, /,'|',divide)                \
                                                                                    \
BINARY_OPERATOR(vectorType, vectorType, tensorType, /,'|',divide)                   \
BINARY_TYPE_OPERATOR(vectorType, vectorType, tensorType, /,'|',divide)              \
                                                                                    \
BINARY_OPERATOR(tensorType, tensorType, tensorType, /,'|',divide)                   \
BINARY_TYPE_OPERATOR(tensorType, tensorType, tensorType, /,'|',divide)              \
                                                                                    \
BINARY_OPERATOR(tensorType, tensorType, diagTensorType, /,'|',divide)               \
BINARY_TYPE_OPERATOR(tensorType, tensorType, diagTensorType, /,'|',divide)          \
                                                                                    \
BINARY_OPERATOR(tensorType, diagTensorType, tensorType, /,'|',divide)               \
BINARY_TYPE_OPERATOR(tensorType, diagTensorType, tensorType, /,'|',divide)          \
                                                                                    \
BINARY_OPERATOR(tensorType, sphericalTensorType, tensorType, /,'|',divide)          \
BINARY_TYPE_OPERATOR(tensorType, sphericalTensorType, tensorType, /,'|',divide)     \
                                                                                    \
BINARY_OPERATOR(tensorType, tensorType, sphericalTensorType, /,'|',divide)          \
BINARY_TYPE_OPERATOR(tensorType, tensorType, sphericalTensorType, /,'|',divide)     \
                                                                                    \
BINARY_OPERATOR(tensorType, tensorType, tensorType, +,'+',add)                      \
BINARY_OPERATOR(tensorType, tensorType, tensorType, -,'-',subtract)                 \
                                                                                    \
BINARY_TYPE_OPERATOR(tensorType, tensorType, tensorType, +,'+',add)                 \
BINARY_TYPE_OPERATOR(tensorType, tensorType, tensorType, -,'-',subtract)            \
                                                                                    \
BINARY_OPERATOR(tensorType, diagTensorType, tensorType, +,'+',add)                  \
BINARY_OPERATOR(tensorType, diagTensorType, tensorType, -,'-',subtract)             \
                                                                                    \
BINARY_TYPE_OPERATOR(tensorType, diagTensorType, tensorType, +,'+',add)             \
BINARY_TYPE_OPERATOR(tensorType, diagTensorType, tensorType, -,'-',subtract)        \
                                                                                    \
BINARY_OPERATOR(tensorType, sphericalTensorType, tensorType, +,'+',add)             \
BINARY_OPERATOR(tensorType, sphericalTensorType, tensorType, -,'-',subtract)        \
                                                                                    \
BINARY_TYPE_OPERATOR(tensorType, sphericalTensorType, tensorType, +,'+',add)        \
BINARY_TYPE_OPERATOR(tensorType, sphericalTensorType, tensorType, -,'-',subtract)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

forVectorTensor5Types(TensorN_FieldFunctions)

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef TensorN_FieldFunctions

#include "undefFieldFunctionsM.H"

#define TEMPLATE template<template<class> class PatchField, class GeoMesh>
#include "GeometricFieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define DiagTensorN_FieldFunctions(tensorType, diagTensorType, sphericalTensorType,     \
    vectorType, CmptType, args...)                                                      \
                                                                                        \
UNARY_FUNCTION(diagTensorType, diagTensorType,inv,inv)                                  \
UNARY_FUNCTION(diagTensorType, diagTensorType,diag,diag)                                \
UNARY_FUNCTION(CmptType, diagTensorType,contractScalar,contractScalar)                  \
UNARY_FUNCTION(vectorType, diagTensorType,contractLinear,contractLinear)                \
                                                                                        \
BINARY_OPERATOR(diagTensorType, CmptType, diagTensorType, *,'*',multiply)               \
BINARY_TYPE_OPERATOR(diagTensorType, CmptType, diagTensorType, *,'*',multiply)          \
                                                                                        \
BINARY_OPERATOR(diagTensorType, CmptType, diagTensorType, /,'|',divide)                 \
BINARY_TYPE_OPERATOR(diagTensorType, CmptType, diagTensorType, /,'|',divide)            \
                                                                                        \
BINARY_OPERATOR(vectorType, vectorType, diagTensorType, /,'|',divide)                   \
BINARY_TYPE_OPERATOR(vectorType, vectorType, diagTensorType, /,'|',divide)              \
                                                                                        \
BINARY_OPERATOR(diagTensorType, diagTensorType, diagTensorType, /,'|',divide)           \
BINARY_TYPE_OPERATOR(diagTensorType, diagTensorType, diagTensorType, /,'|',divide)      \
                                                                                        \
BINARY_OPERATOR(diagTensorType, sphericalTensorType, diagTensorType, /,'|',divide)      \
BINARY_TYPE_OPERATOR(diagTensorType, sphericalTensorType, diagTensorType, /,'|',divide) \
                                                                                        \
BINARY_OPERATOR(diagTensorType, diagTensorType, sphericalTensorType, /,'|',divide)      \
BINARY_TYPE_OPERATOR(diagTensorType, diagTensorType, sphericalTensorType, /,'|',divide) \
                                                                                        \
BINARY_OPERATOR(diagTensorType, diagTensorType, diagTensorType, +,'+',add)              \
BINARY_OPERATOR(diagTensorType, diagTensorType, diagTensorType, -,'-',subtract)         \
                                                                                        \
BINARY_TYPE_OPERATOR(diagTensorType, diagTensorType, diagTensorType, +,'+', add)        \
BINARY_TYPE_OPERATOR(diagTensorType, diagTensorType, diagTensorType, -,'-', subtract)   \
                                                                                        \
BINARY_OPERATOR(diagTensorType, sphericalTensorType, diagTensorType, +,'+', add)        \
BINARY_OPERATOR(diagTensorType, sphericalTensorType, diagTensorType, -,'-', subtract)   \
                                                                                        \
BINARY_TYPE_OPERATOR(diagTensorType, sphericalTensorType, diagTensorType, +,'+', add)   \
BINARY_TYPE_OPERATOR(diagTensorType, sphericalTensorType, diagTensorType, -,'-', subtract)  \
                                                                                        \
BINARY_OPERATOR(diagTensorType, diagTensorType, sphericalTensorType, +,'+', add)        \
BINARY_OPERATOR(diagTensorType, diagTensorType, sphericalTensorType, -,'-', subtract)   \
                                                                                        \
BINARY_TYPE_OPERATOR(diagTensorType, diagTensorType, sphericalTensorType, +,'+', add)   \
BINARY_TYPE_OPERATOR(diagTensorType, diagTensorType, sphericalTensorType, -,'+', subtract)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

forVectorTensor5Types(DiagTensorN_FieldFunctions)

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef DiagTensorN_FieldFunctions

#include "undefFieldFunctionsM.H"

#define TEMPLATE template<template<class> class PatchField, class GeoMesh>
#include "GeometricFieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define SphericalTensorN_FieldFunctions(tensorType, diagTensorType,                 \
    sphericalTensorType, vectorType, CmptType, args...)                             \
                                                                                    \
UNARY_FUNCTION(sphericalTensorType, sphericalTensorType,inv,inv)                    \
UNARY_FUNCTION(sphericalTensorType, sphericalTensorType,diag,diag)                  \
                                                                                    \
BINARY_OPERATOR(sphericalTensorType, CmptType, sphericalTensorType, *,'*',multiply)                   \
BINARY_TYPE_OPERATOR(sphericalTensorType, CmptType, sphericalTensorType, *,'*',multiply)              \
                                                                                    \
BINARY_OPERATOR(sphericalTensorType, CmptType, sphericalTensorType, /,'|',divide)   \
BINARY_TYPE_OPERATOR(sphericalTensorType, CmptType, sphericalTensorType, /,'|',divide) \
                                                                                    \
BINARY_OPERATOR(vectorType, vectorType, sphericalTensorType, /,'|',divide)          \
BINARY_TYPE_OPERATOR(vectorType, vectorType, sphericalTensorType, /,'|',divide)     \
                                                                                    \
BINARY_OPERATOR(sphericalTensorType, sphericalTensorType, sphericalTensorType, /,'|',divide)          \
BINARY_TYPE_OPERATOR(sphericalTensorType, sphericalTensorType, sphericalTensorType, /,'|',divide)     \
                                                                                    \
BINARY_OPERATOR(sphericalTensorType, sphericalTensorType, sphericalTensorType, +,'+',add)             \
BINARY_OPERATOR(sphericalTensorType, sphericalTensorType, sphericalTensorType, -,'-',subtract)        \
                                                                                    \
BINARY_TYPE_OPERATOR(sphericalTensorType, sphericalTensorType, sphericalTensorType, +,'+',add)        \
BINARY_TYPE_OPERATOR(sphericalTensorType, sphericalTensorType, sphericalTensorType, -,'-',subtract)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

forVectorTensor5Types(SphericalTensorN_FieldFunctions)

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef SphericalTensorN_FieldFunctions

#include "undefFieldFunctionsM.H"

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#define doMakeVolFields(type, Type, args...)                    \
    defineTemplateTypeNameAndDebug(                             \
        vol##Type##Field::DimensionedInternalField, 0);         \
                                                                \
    defineTemplateTypeNameAndDebug(vol##Type##Field, 0);

forVector5Type(doMakeVolFields)

forTensor5Type(doMakeVolFields)

forDiagTensor5Type(doMakeVolFields)

forSphericalTensor5Type(doMakeVolFields)

#undef doMakeVolFields

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam



#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#define doMakeSurfaceFields(type, Type, args...)                \
    defineTemplateTypeNameAndDebug(                             \
        surface##Type##Field::DimensionedInternalField, 0);     \
                                                                \
    defineTemplateTypeNameAndDebug(surface##Type##Field, 0);

forVector5Type(doMakeSurfaceFields)

forTensor5Type(doMakeSurfaceFields)

forDiagTensor5Type(doMakeSurfaceFields)

forSphericalTensor5Type(doMakeSurfaceFields)

#undef doMakeSurfaceFields

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "fvPatchField.H"
#include "fvPatchFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#define makeFvPatchField(fvPatchTypeField)                                    \
                                                                              \
defineNamedTemplateTypeNameAndDebug(fvPatchTypeField, 0);                     \
template<>                                                                    \
debug::debugSwitch                                                            \
fvPatchTypeField::disallowGenericFvPatchField                                 \
(                                                                             \
    "disallowGenericFvPatchField",                                            \
    0                                                                         \
);                                                                            \
defineTemplateRunTimeSelectionTable(fvPatchTypeField, patch);                 \
defineTemplateRunTimeSelectionTable(fvPatchTypeField, patchMapper);           \
defineTemplateRunTimeSelectionTable(fvPatchTypeField, dictionary);


#define doMakeFvPatchField(type, Type, args...)                               \
    makeFvPatchField(fvPatch##Type##Field)

forVector5Type(doMakeFvPatchField)

forTensor5Type(doMakeFvPatchField)

forDiagTensor5Type(doMakeFvPatchField)

forSphericalTensor5Type(doMakeFvPatchField)

#undef doMakeFvPatchField
#undef makeFvPatchField

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "calculatedFvPatchField.H"
#include "volMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#define doMakePatchTypeField(type, Type, args...)                             \
                                                                              \
makeTemplatePatchTypeField                                                    \
(                                                                             \
    fvPatch##Type##Field,                                                     \
    calculatedFvPatch##Type##Field                                            \
);

forVector5Type(doMakePatchTypeField)

forTensor5Type(doMakePatchTypeField)

forDiagTensor5Type(doMakePatchTypeField)

forSphericalTensor5Type(doMakePatchTypeField)

#undef doMakePatchTypeField

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "coupledFvPatchVectorNFieldsFwd.H"
#include "coupledFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#define doMakePatchTypeField(type, Type, args...)           \
    makePatchTypeFieldTypeName(coupledFvPatch##Type##Field);

forVector5Type(doMakePatchTypeField)

forTensor5Type(doMakePatchTypeField)

forDiagTensor5Type(doMakePatchTypeField)

forSphericalTensor5Type(doMakePatchTypeField)


#undef doMakePatchTypeField

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "volVectorNFields.H"
#include "processorFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#define doMakePatchTypeField(type, Type, args...)                             \
                                                                              \
makeTemplatePatchTypeField                                                    \
(                                                                             \
    fvPatch##Type##Field,                                                     \
    processorFvPatch##Type##Field                                             \
);

forVector5Type(doMakePatchTypeField)

#undef doMakePatchTypeField

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "fvsPatchField.H"
#include "fvsPatchFieldsFwd.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#define makeFvsPatchField(fvsPatchTypeField)                                  \
                                                                              \
defineNamedTemplateTypeNameAndDebug(fvsPatchTypeField, 0);                    \
template<>                                                                    \
Foam::debug::debugSwitch                                                      \
fvsPatchTypeField::disallowDefaultFvsPatchField                               \
(                                                                             \
    "disallowDefaultFvsPatchField",                                           \
    0                                                                         \
);                                                                            \
defineTemplateRunTimeSelectionTable(fvsPatchTypeField, patch);                \
defineTemplateRunTimeSelectionTable(fvsPatchTypeField, patchMapper);          \
defineTemplateRunTimeSelectionTable(fvsPatchTypeField, dictionary);

#define doMakeFvsPatchField(type, Type, args...)    \
    makeFvsPatchField(fvsPatch##Type##Field)

forVector5Type(doMakeFvsPatchField)

forTensor5Type(doMakeFvsPatchField)

forDiagTensor5Type(doMakeFvsPatchField)

forSphericalTensor5Type(doMakeFvsPatchField)

#undef makeFvsPatchField
#undef doMakeFvsPatchField

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "calculatedFvsPatchField.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#define doMakeFvsPatchTypeField(type, Type, args...)                          \
    makeFvsPatchTypeField(fvsPatch##Type##Field, calculatedFvsPatch##Type##Field);

forVector5Type(doMakeFvsPatchTypeField)

forTensor5Type(doMakeFvsPatchTypeField)

forDiagTensor5Type(doMakeFvsPatchTypeField)

forSphericalTensor5Type(doMakeFvsPatchTypeField)

#undef doMakeFvsPatchTypeField

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "BlockLduInterfaceField.H"
#include "GGIBlockLduInterfaceField.H"
#include "MixingPlaneBlockLduInterfaceField.H"
#include "OverlapGGIBlockLduInterfaceField.H"
#include "ProcessorBlockLduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#define makeTemplateTypeNameAndDebug(type, Type, args...)                     \
                                                                              \
defineTemplateTypeNameAndDebug(BlockLduInterfaceField<type>, 0);              \
                                                                              \
defineTemplateTypeNameAndDebug(GGIBlockLduInterfaceField<type>, 0);           \
defineTemplateTypeNameAndDebug(MixingPlaneBlockLduInterfaceField<type>, 0);   \
defineTemplateTypeNameAndDebug(OverlapGGIBlockLduInterfaceField<type>, 0);    \
defineTemplateTypeNameAndDebug(ProcessorBlockLduInterfaceField<type>, 0);

forVector5Type(makeTemplateTypeNameAndDebug);

#undef makeTemplateTypeNameAndDebug

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "coeffFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#define makeNamedTemplateTypeNameAndDebug(type, Type, args...)                \
    defineNamedTemplateTypeNameAndDebug(block##Type##Matrix, 0);

forVector5Type(makeNamedTemplateTypeNameAndDebug);

#undef makeNamedTemplateTypeNameAndDebug

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineNamedTemplateTypeNameAndDebug(blockAmgVector5Level, 0);

defineNamedTemplateTypeNameAndDebug(coarseblockAmgVector5Level, 0);

defineNamedTemplateTypeNameAndDebug(fineblockAmgVector5Level, 0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineNamedTemplateTypeNameAndDebug(blockAMGVector5Cycle, 0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "ProcessorBlockSAMGInterfaceField.H"
#include "GGIBlockSAMGInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#define makeTemplateTypeNameAndDebug(type, Type, args...)                     \
                                                                              \
typedef BlockSAMGInterfaceField<type > block##Type##SAMGInterfaceField;       \
defineNamedTemplateTypeNameAndDebug(block##Type##SAMGInterfaceField, 0);      \
defineTemplateRunTimeSelectionTable(block##Type##SAMGInterfaceField, lduInterface); \
                                                                              \
typedef ProcessorBlockSAMGInterfaceField<type > block##Type##ProcessorSAMGInterfaceField;  \
makeBlockSAMGInterfaceField(block##Type##SAMGInterfaceField, block##Type##ProcessorSAMGInterfaceField); \
                                                                              \
typedef GGIBlockSAMGInterfaceField<type > block##Type##GGISAMGInterfaceField; \
makeBlockSAMGInterfaceField(block##Type##SAMGInterfaceField, block##Type##GGISAMGInterfaceField); \

forVector5Type(makeTemplateTypeNameAndDebug);

#undef makeTemplateTypeNameAndDebug

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "ProcessorBlockAMGInterfaceField.H"
#include "GGIBlockAMGInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#define makeTemplateTypeNameAndDebug(type, Type, args...)                     \
                                                                              \
typedef BlockAMGInterfaceField<type > block##Type##AMGInterfaceField;         \
defineNamedTemplateTypeNameAndDebug(block##Type##AMGInterfaceField, 0);       \
defineTemplateRunTimeSelectionTable(block##Type##AMGInterfaceField, lduInterface); \
                                                                              \
typedef ProcessorBlockAMGInterfaceField<type > block##Type##ProcessorAMGInterfaceField;  \
makeBlockAMGInterfaceField(block##Type##AMGInterfaceField, block##Type##ProcessorAMGInterfaceField); \
typedef GGIBlockAMGInterfaceField<type > block##Type##GGIAMGInterfaceField;  \
makeBlockAMGInterfaceField(block##Type##AMGInterfaceField, block##Type##GGIAMGInterfaceField); \

forVector5Type(makeTemplateTypeNameAndDebug);

#undef makeTemplateTypeNameAndDebug

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "blockLduMatrices.H"

#include "blockLduPrecons.H"
#include "BlockNoPrecon.H"
#include "blockDiagonalPrecons.H"
#include "blockGaussSeidelPrecons.H"
#include "BlockCholeskyPrecon.H"
#include "blockILUC0Precons.H"
#include "BlockILUCpPrecon.H"

#include "blockLduSmoothers.H"
#include "blockGaussSeidelSmoothers.H"
#include "BlockILUSmoother.H"
#include "blockILUC0Smoothers.H"
#include "BlockILUCpSmoother.H"

#include "blockLduSolvers.H"
#include "BlockDiagonalSolver.H"
#include "BlockBiCGStabSolver.H"
#include "BlockCGSolver.H"
#include "BlockGaussSeidelSolver.H"
#include "BlockILUSolver.H"
#include "BlockGMRESSolver.H"

#include "blockAMGSolvers.H"
#include "blockAMGPrecons.H"
#include "blockMatrixCoarsenings.H"
#include "blockMatrixClusterings.H"
#include "blockMatrixSelections.H"
#include "blockCoeffNorms.H"
#include "blockCoeffTwoNorms.H"
#include "blockCoeffMaxNorms.H"
#include "blockCoeffComponentNorms.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#define makeSolver(type, Type, args...)                                       \
/* Preconditioners */                                                         \
typedef BlockLduPrecon<type > block##Type##Precon;                            \
defineNamedTemplateTypeNameAndDebug(block##Type##Precon, 0);                  \
defineTemplateRunTimeSelectionTable(block##Type##Precon, dictionary);         \
                                                                              \
typedef BlockNoPrecon<type > block##Type##NoPrecon;                           \
makeBlockPrecon(block##Type##Precon, block##Type##NoPrecon);                  \
                                                                              \
typedef BlockDiagonalPrecon<type > block##Type##DiagonalPrecon;               \
makeBlockPrecon(block##Type##Precon, block##Type##DiagonalPrecon);            \
                                                                              \
typedef BlockGaussSeidelPrecon<type > block##Type##GaussSeidelPrecon;         \
makeBlockPrecon(block##Type##Precon, block##Type##GaussSeidelPrecon);         \
                                                                              \
typedef BlockCholeskyPrecon<type > block##Type##CholeskyPrecon;               \
makeBlockPrecon(block##Type##Precon, block##Type##CholeskyPrecon);            \
                                                                              \
typedef BlockILUC0Precon<type > block##Type##ILUC0Precon;                     \
makeBlockPrecon(block##Type##Precon, block##Type##ILUC0Precon);               \
                                                                              \
typedef BlockILUCpPrecon<type > block##Type##ILUCpPrecon;                     \
makeBlockPrecon(block##Type##Precon, block##Type##ILUCpPrecon);               \
                                                                              \
/* Smoothers */                                                               \
typedef BlockLduSmoother<type > block##Type##Smoother;                        \
defineNamedTemplateTypeNameAndDebug(block##Type##Smoother, 0);                \
defineTemplateRunTimeSelectionTable(block##Type##Smoother, dictionary);       \
                                                                              \
typedef BlockGaussSeidelSmoother<type > block##Type##GaussSeidelSmoother;     \
makeBlockSmoother(block##Type##Smoother, block##Type##GaussSeidelSmoother);   \
                                                                              \
typedef BlockILUSmoother<type > block##Type##ILUSmoother;                     \
makeBlockSmoother(block##Type##Smoother, block##Type##ILUSmoother);           \
                                                                              \
typedef BlockILUC0Smoother<type > block##Type##ILUC0Smoother;                 \
makeBlockSmoother(block##Type##Smoother, block##Type##ILUC0Smoother);         \
                                                                              \
typedef BlockILUCpSmoother<type > block##Type##ILUCpSmoother;                 \
makeBlockSmoother(block##Type##Smoother, block##Type##ILUCpSmoother);         \
                                                                              \
                                                                              \
/* Solvers */                                                                 \
typedef BlockLduSolver<type > block##Type##Solver;                            \
defineNamedTemplateTypeNameAndDebug(block##Type##Solver, 0);                  \
defineTemplateRunTimeSelectionTable                                           \
(                                                                             \
    block##Type##Solver,                                                      \
    symMatrix                                                                 \
);                                                                            \
                                                                              \
defineTemplateRunTimeSelectionTable                                           \
(                                                                             \
    block##Type##Solver,                                                      \
    asymMatrix                                                                \
);                                                                            \
                                                                              \
typedef BlockDiagonalSolver<type > block##Type##DiagonalSolver;               \
defineNamedTemplateTypeNameAndDebug(block##Type##DiagonalSolver, 0);          \
                                                                              \
typedef BlockBiCGStabSolver<type > block##Type##BiCGStabSolver;               \
makeBlockSolverTypeName(block##Type##BiCGStabSolver);                         \
addSolverToBlockMatrix(Type, block##Type##BiCGStabSolver, symMatrix);         \
addSolverToBlockMatrix(Type, block##Type##BiCGStabSolver, asymMatrix);        \
                                                                              \
typedef BlockCGSolver<type > block##Type##CGSolver;                           \
makeBlockSolverTypeName(block##Type##CGSolver);                               \
addSolverToBlockMatrix(Type, block##Type##CGSolver, symMatrix);               \
                                                                              \
typedef BlockGaussSeidelSolver<type > block##Type##GaussSeidelSolver;         \
makeBlockSolverTypeName(block##Type##GaussSeidelSolver);                      \
addSolverToBlockMatrix(Type, block##Type##GaussSeidelSolver, symMatrix);      \
addSolverToBlockMatrix(Type, block##Type##GaussSeidelSolver, asymMatrix);     \
                                                                              \
typedef BlockILUSolver<type > block##Type##ILUSolver;                         \
makeBlockSolverTypeName(block##Type##ILUSolver);                              \
addSolverToBlockMatrix(Type, block##Type##ILUSolver, symMatrix);              \
addSolverToBlockMatrix(Type, block##Type##ILUSolver, asymMatrix);             \
                                                                              \
typedef BlockGMRESSolver<type > block##Type##GMRESSolver;                     \
makeBlockSolverTypeName(block##Type##GMRESSolver);                            \
addSolverToBlockMatrix(Type, block##Type##GMRESSolver, symMatrix);            \
addSolverToBlockMatrix(Type, block##Type##GMRESSolver, asymMatrix);           \
                                                                              \
typedef BlockMatrixCoarsening<type > block##Type##MatrixCoarsening;           \
defineNamedTemplateTypeNameAndDebug(block##Type##MatrixCoarsening, 0);        \
defineTemplateRunTimeSelectionTable(block##Type##MatrixCoarsening, matrix);   \
                                                                              \
typedef BlockMatrixClustering<type > block##Type##MatrixClustering;           \
makeBlockMatrixCoarsening(block##Type##MatrixCoarsening, block##Type##MatrixClustering); \
                                                                              \
typedef BlockMatrixSelection<type > block##Type##MatrixSelection;             \
makeBlockMatrixCoarsening(block##Type##MatrixCoarsening, block##Type##MatrixSelection); \
                                                                              \
typedef BlockCoeffNorm<type > block##Type##CoeffNorm;                         \
defineNamedTemplateTypeNameAndDebug(block##Type##CoeffNorm, 0);               \
defineTemplateRunTimeSelectionTable(block##Type##CoeffNorm, dictionary);      \
                                                                              \
typedef BlockCoeffTwoNorm<type > block##Type##CoeffTwoNorm;                   \
makeBlockCoeffNorm(block##Type##CoeffNorm, block##Type##CoeffTwoNorm);        \
                                                                              \
typedef BlockCoeffComponentNorm<type > block##Type##CoeffComponentNorm;       \
makeBlockCoeffNorm(block##Type##CoeffNorm, block##Type##CoeffComponentNorm);  \
                                                                              \
typedef BlockCoeffMaxNorm<type > block##Type##CoeffMaxNorm;                   \
makeBlockCoeffNorm(block##Type##CoeffNorm, block##Type##CoeffMaxNorm);        \
                                                                              \
typedef BlockAMGSolver<type > block##Type##AmgSolver;                         \
makeBlockSolverTypeName(block##Type##AmgSolver);                              \
addSolverToBlockMatrix(Type, block##Type##AmgSolver, symMatrix);              \
addSolverToBlockMatrix(Type, block##Type##AmgSolver, asymMatrix);             \
                                                                              \
typedef BlockAMGPrecon<type > block##Type##AmgPrecon;                         \
makeBlockPrecon(block##Type##Precon, block##Type##AmgPrecon);                 \

forVector5Type(makeSolver)

#undef makeSolver

} // End namespace Foam

#undef forVector5Type
#undef forTensor5Type
#undef forDiagTensor5Type
#undef forSphericalTensor5Type
#undef forVectorTensor5Types

// ************************************************************************* //
