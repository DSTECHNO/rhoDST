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
    rhoDST::jacobian

Description
    Stores and inserts ownerCoeffs and neighbourCoeffs. 

SourceFiles
    jacobian.C
    
\*---------------------------------------------------------------------------------------*/

#ifndef jacobian_H
#define jacobian_H

#include "IOdictionary.H"
#include "volFields.H"
#include "fvCFD.H"
#include "Vector5typedefs.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class jacobian Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class jacobian
{
    private:

    protected:
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> J_;
    
    public:

            //- Constructors
            jacobian(const fvMesh& mesh);
    
        //- Destructor
            ~jacobian();

        //- Member Functions
            GeometricField<Type, fvsPatchField, surfaceMesh>& J();
            typename Foam::GeometricField<Type, fvsPatchField, surfaceMesh>::GeometricBoundaryField& boundaryField();

            template<class fType>
            void insertScalarEqn
            (
                tmp<GeometricField<fType, fvsPatchField, surfaceMesh>> tJi,
                int row,
                int col    
            );
        
            template<class fType>
            void insertScalarEqn
            (
                const GeometricField<fType, fvsPatchField, surfaceMesh>& Ji,
                int row,
                int col    
            );
        
            template<class fType>
            void insertVectorEqn
            (
                tmp<GeometricField<fType, fvsPatchField, surfaceMesh>> tJi,
                int row,
                int col    
            );
        
            template<class fType>
            void insertVectorEqn
            (
                const GeometricField<fType, fvsPatchField, surfaceMesh>& Ji,
                int row,
                int col    
            );
    
           //- Operators
            tmp<Field<Type>> operator*(const Field<scalar>& Sf) const;
            void operator*=(const surfaceScalarField& w);
    
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
}

#ifdef NoRepository
#   include "jacobian.C"
#endif

#endif
