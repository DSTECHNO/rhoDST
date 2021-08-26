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

#include "fvCFD.H"
#include "fvMesh.H"
#include "fluxDSTIntegration.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvmDST
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<fvBlockMatrixDST<Type>> fluxDST
(
    baseFluxDST& rhoDSTFlux,
    geometricFieldDST<Type>& psi
)
{
    return fvDST::fluxDSTIntegration<Type>::New
            (
                psi.dictDST()
            )().integrate
            (
                rhoDSTFlux,
                psi
            );
}


template<class Type>
tmp<fvBlockMatrixDST<Type>> ddt
(
    geometricFieldDST<Type>& psi
)
{
        const fvMesh& mesh = psi.mesh();        

        tmp<fvBlockMatrixDST<Type>> tEqn(new fvBlockMatrixDST<Type>(psi));
        fvBlockMatrixDST<Type>& Eqn = tEqn();

        typename CoeffField<Type>::squareTypeField& d = Eqn.diag().asSquare();

        typename CoeffField<Type>::squareTypeField& u = Eqn.upper().asSquare();

        typename CoeffField<Type>::squareTypeField& l = Eqn.lower().asSquare();
        Field<Type>& s = Eqn.source();

        // The matrix diagonal coeffs are calculated using the ddtScheme for rho variable given in
        // fvSchemes dict. If rho, rhoU and rhoE have the same ddtScheme then, the diagonals would
        // be the same for each variable. 
        tmp<fvScalarMatrix> rhoEqn(fvm::ddt(psi.rho()));

        d = DST::I<Type>()*rhoEqn->diag();
        DST::insert(s, rhoEqn->source()-rhoEqn->diag()*psi.rho().internalField(), 0);
        rhoEqn.clear();
        
        tmp<fvVectorMatrix> rhoUEqn(fvm::ddt(psi.rhoU()));
        DST::insert(s, rhoUEqn->source()-rhoUEqn->diag()*psi.rhoU().internalField(), 1);
        rhoUEqn.clear();
        
        tmp<fvScalarMatrix> rhoEEqn(fvm::ddt(psi.rhoE()));
        DST::insert(s, rhoEEqn->source()-rhoEEqn->diag()*psi.rhoE().internalField(), 4);
        rhoEEqn.clear();    
        
        
        forAll(psi.boundaryField(), patchI)
        {
            if (psi.boundaryField()[patchI].patch().coupled())
            {
                typename CoeffField<Type>::squareTypeField& pCoupleUpper = Eqn.coupleUpper()[patchI].asSquare();
                typename CoeffField<Type>::squareTypeField& pCoupleLower = Eqn.coupleLower()[patchI].asSquare();
            }
        }    
    
        return  tEqn;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvm

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
