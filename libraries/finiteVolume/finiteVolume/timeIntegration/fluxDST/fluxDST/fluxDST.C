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

#include "fvCFD.H"
#include "fluxDST.H"
#include "addToRunTimeSelectionTable.H"
#include "fvBlockMatrix.H"
#include "BlockLduSystem.H"
#include "DST.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fvDST
{

defineNamedTemplateTypeNameAndDebug(fluxDST<vector5>, 0);
addTemplatedToRunTimeSelectionTable(fluxDSTIntegration, fluxDST, vector5, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
fluxDST<Type>::fluxDST()
{} // End of constructor

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
fluxDST<Type>::~fluxDST()
{} //End of destructor

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
template<class Type>
tmp<fvBlockMatrixDST<Type>> fluxDST<Type>::integrate
(
    baseFluxDST& rhoDSTFlux,
    geometricFieldDST<Type>& W
) const 
{
    JacobianDST<typename outerProduct< Type, Type >::type>& dRdW = rhoDSTFlux.Jacobian();
    const fvMesh& mesh = W.mesh();

    tmp<fvBlockMatrixDST<Type>> fEqn(new fvBlockMatrixDST<Type>(W));
    fvBlockMatrixDST<Type>& Eqn = fEqn();
    
    typename CoeffField<Type>::squareTypeField& d = Eqn.diag().asSquare();
    typename CoeffField<Type>::squareTypeField& u = Eqn.upper().asSquare();
    typename CoeffField<Type>::squareTypeField& l = Eqn.lower().asSquare();
    Field<Type>& s = Eqn.source();
    
    
    // beta = 1 by default, used for dual time integration
    scalar beta = W.dictDST().beta();
    
    const vectorField& SfIn = mesh.Sf();
    const scalarField& magSfIn = mesh.magSf();    

    l = -beta*(dRdW.ownerCoeffs()*magSfIn + dRdW.spectInternal());
    u = beta*(dRdW.neighbourCoeffs()*magSfIn - dRdW.spectInternal());
    Eqn.negSumDiag();    
    
    forAll(W.boundaryField(), patchI)
    {
        const fvPatch& patch = W.boundaryField()[patchI].patch();
        const scalarField& magSf = patch.magSf();
        const unallocLabelList& fc = patch.faceCells();
        
        Field<typename outerProduct< Type, Type >::type> pd(beta*(dRdW.patchInternalCoeffs()[patchI] + dRdW.spectPatch()[patchI]));
        forAll(pd, faceI)
        {
            d[fc[faceI]] += pd[faceI];
        }
        
        if (patch.coupled())
        {
            typename CoeffField<Type>::squareTypeField& pCoupleUpper = Eqn.coupleUpper()[patchI].asSquare();
            typename CoeffField<Type>::squareTypeField& pCoupleLower = Eqn.coupleLower()[patchI].asSquare();

            pCoupleLower = beta*(dRdW.patchInternalCoeffs()[patchI] + dRdW.spectPatch()[patchI]);
            pCoupleUpper = beta*(-dRdW.patchNeighbourCoeffs()[patchI] + dRdW.spectPatch()[patchI]);
        }
    } 
    
    s -= rhoDSTFlux.R()*mesh.V();
        
    return fEqn;
}

}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
}

