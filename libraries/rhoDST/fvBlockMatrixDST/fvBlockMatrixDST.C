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
template <class Type>
fvBlockMatrixDST<Type>::fvBlockMatrixDST
(
    const BlockLduSystem<Type, Type>& bM,
    geometricFieldDST<Type>& bP
)
:
    BlockLduSystem<Type, Type>(bM),
    psi_(bP)
{}

template <class Type>
fvBlockMatrixDST<Type>::fvBlockMatrixDST
(
    geometricFieldDST<Type>& bP
)
:
    BlockLduSystem<Type, Type>(bP.mesh()),
    psi_(bP)
{}

template <class Type>
fvBlockMatrixDST<Type>::fvBlockMatrixDST
(
    const fvBlockMatrixDST<Type>& fvBM
)
:
    BlockLduSystem<Type, Type>(fvBM),
    psi_(fvBM.psi())
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
template <class Type>
fvBlockMatrixDST<Type>::~fvBlockMatrixDST()
{}


// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

template<class Type>
geometricFieldDST<Type>& fvBlockMatrixDST<Type>::psi() const
{
    return psi_;
}

template<class Type>
geometricFieldDST<Type>& fvBlockMatrixDST<Type>::psi()
{
    return psi_;
}

template<class Type>
Foam::BlockSolverPerformance<Type> fvBlockMatrixDST<Type>::solve
(
    bool print
)
{
    
    this->interfaces() = psi_.boundaryField().blockInterfaces();
    word psiName = psi_.dictDST().template lookupOrDefault<word>("dependentVariable","W");
    if(psiName=="W")
    {
        dictionary solverControls(psi_.mesh().solutionDict().solverDict(psi_.name()));


        Field<Type> sourceW(psi_.internalField().size(), Type::zero);
        this->Amul(sourceW, psi_.internalField());
        sourceW += this->source();

      
        BlockSolverPerformance<Type> solverPerf =
            BlockLduSolver<Type>::New
            (
                psi_.name(),
                *this,
                solverControls
            )->solve(psi_.internalField(), sourceW);


        if(print)
            solverPerf.print();


        psi_.mesh().solutionDict().setSolverPerformance(psi_.name(), solverPerf);

        return solverPerf;
    }
    else
    {
        dictionary solverControls(psi_.mesh().solutionDict().solverDict(psi_.delta().name()));
      
        BlockSolverPerformance<Type> solverPerf =
            BlockLduSolver<Type>::New
            (
                psi_.delta().name(),
                *this,
                solverControls
            )->solve(psi_.delta().internalField(), this->source());
        psi_++;

        if(print)
            solverPerf.print();


        psi_.mesh().solutionDict().setSolverPerformance(psi_.name(), solverPerf);

        return solverPerf;
    }
}

template<class Type>
tmp<Field<Type>> fvBlockMatrixDST<Type>::residualDST()
{
    tmp<Field<Type>> tR(new Field<Type>(this->source()/psi_.mesh().V()));
    return tR;
}

// * * * * * * * * * * * * * * * * Operators  * * * * * * * * * * * * * * * //

template<class Type>
fvBlockMatrixDST<Type> fvBlockMatrixDST<Type>::operator+
(
    const fvBlockMatrixDST<Type>& mB
)
{
    BlockLduSystem<Type, Type>::operator+=(mB);
    return (*this);
}


template<class Type>
fvBlockMatrixDST<Type> fvBlockMatrixDST<Type>::operator+
(
    const vector5Field& sB
)
{
    Field<Type>& s(this->source());
    s -= sB;
    
    return (*this);
}


template<class Type>
void fvBlockMatrixDST<Type>::operator+=
(
    const fvBlockMatrixDST<Type>& mB
)
{
    BlockLduSystem<Type, Type>::operator+=(mB);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template <class Type>
tmp<fvBlockMatrixDST<Type>> operator+
(
    const tmp<fvBlockMatrixDST<Type>>& a, 
    const tmp<fvBlockMatrixDST<Type>>& b
)
{

    tmp<fvBlockMatrixDST<Type>> c(a);
    c() += b();
    b.clear();
    return c;
}

} // End namespace Foam
