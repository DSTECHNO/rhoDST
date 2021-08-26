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

\*---------------------------------------------------------------------------------------*/

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
jacobian<Type>::jacobian
(
    const fvMesh& mesh
)
:
    J_
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "J_",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>("J_",dimless,Type::zero)
        )
    )
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
jacobian<Type>::~jacobian()
{}

// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

template<class Type>
GeometricField<Type, fvsPatchField, surfaceMesh>& jacobian<Type>::J()
{
    return J_();
}

template<class Type>
typename Foam::GeometricField<Type, fvsPatchField, surfaceMesh>::GeometricBoundaryField& jacobian<Type>::boundaryField()
{
    return J_().boundaryField();
}

template<class Type>
template<class fType>
void jacobian<Type>::insertScalarEqn
(
    const GeometricField<fType, fvsPatchField, surfaceMesh>& Ji,
    int row,
    int col
)
{
    int nCols = pow(pTraits<Type>::nComponents,0.5);
    int nComponents = pTraits<fType>::nComponents;
    int dir = nCols*row + col;
    for(int i = 0; i<nComponents; i++)
    {
        J_().replace(dir, Ji.component(i));
        dir++;
    }
}

template<class Type>
template<class fType>
void jacobian<Type>::insertVectorEqn
(
    const GeometricField<fType, fvsPatchField, surfaceMesh>& Ji,
    int row,
    int col
)
{
    int nCols = pow(pTraits<Type>::nComponents,0.5);
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
            J_().replace(dir, Ji.component(j));
            dir += nCols;
        }
        dir -= nCols*rowF -1;
    }
}

template<class Type>
template<class fType>
void jacobian<Type>::insertScalarEqn
(
    tmp<GeometricField<fType, fvsPatchField, surfaceMesh>> tJi,
    int row,
    int col
)
{
    insertScalarEqn(tJi(), row, col);
}

template<class Type>
template<class fType>
void jacobian<Type>::insertVectorEqn
(
    tmp<GeometricField<fType, fvsPatchField, surfaceMesh>> tJi,
    int row,
    int col
)
{
    insertVectorEqn(tJi(), row, col);
}

// * * * * * * * * * * * * * * * * Operators  * * * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type>> jacobian<Type>::operator*(const Field<scalar>& Sf) const
{
    tmp<Field<Type>> fm(new Field<Type>(J_().internalField()*Sf));
    return fm;
}

template<class Type>
void jacobian<Type>::operator*=(const surfaceScalarField& w)
{
    J_() *= w;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam






