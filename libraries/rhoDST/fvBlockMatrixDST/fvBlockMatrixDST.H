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
    rhoDST::fvBlockMatrixDST

Description
    Calculate fluxes using HLL method

SourceFiles
    fvBlockMatrixDST.C

\*---------------------------------------------------------------------------------------*/

#ifndef fvBlockMatrixDST_H
#define fvBlockMatrixDST_H

#include "IOdictionary.H"
#include "volFields.H"
#include "fvCFD.H"
#include "Vector5typedefs.H"
#include "fvBlockMatrix.H"
#include "basicPsiThermo.H"
#include "turbulenceModel.H"
#include "baseFluxDST.H"
#include "JacobianDST.H"

namespace Foam
{
template <class Type>
class geometricFieldDST;
/*---------------------------------------------------------------------------*\
                           Class fvBlockMatrixDST Declaration
\*---------------------------------------------------------------------------*/
template <class Type>
class fvBlockMatrixDST 
:
    public BlockLduSystem<Type, Type>
{            
    public:

        geometricFieldDST<Type>& psi_;
        Type residualDSTV;
    
        //- Constructors
            fvBlockMatrixDST
            (
                const BlockLduSystem<Type, Type>& bM,
                geometricFieldDST<Type>& bP
            );

            fvBlockMatrixDST
            (
                geometricFieldDST<Type>& bP
            );

            fvBlockMatrixDST
            (
                const fvBlockMatrixDST<Type>& fvBM
            );    
    
        //- Destructor
            ~fvBlockMatrixDST();        
    
       //- Member Functions
            geometricFieldDST<Type>& psi() const;
            geometricFieldDST<Type>& psi();
            tmp<Field<Type>> residualDST();

            Foam::BlockSolverPerformance<Type> solve
            (
                bool print
            );
        
        
            //- Operators
            fvBlockMatrixDST<Type> operator+
            (
                const fvBlockMatrixDST<Type>& mB
            );
            void operator+=
            (
                const fvBlockMatrixDST<Type>& mB
            );
            fvBlockMatrixDST<Type> operator+
            (
                const vector5Field& sB
            );
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template <class Type>
tmp<fvBlockMatrixDST<Type>> operator+
(
    const tmp<fvBlockMatrixDST<Type>>& a, 
    const tmp<fvBlockMatrixDST<Type>>& b
);

} // End namespace Foam

#ifdef NoRepository
#   include "fvBlockMatrixDST.C"
#endif

#endif
