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

#include "farFieldTemperatureDSTFvPatchScalarField.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "basicPsiThermo.H"
#include "DST.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

farFieldTemperatureDSTFvPatchScalarField::farFieldTemperatureDSTFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    UInf_(vector::zero),
    pInf_(0),
    TInf_(300)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


farFieldTemperatureDSTFvPatchScalarField::farFieldTemperatureDSTFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    UInf_(dict.lookup("UInf")),
    pInf_(readScalar(dict.lookup("pInf"))),
    TInf_(readScalar(dict.lookup("TInf")))
{
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchScalarField::operator=(this->refValue());
    }

    this->refValue() = *this;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


farFieldTemperatureDSTFvPatchScalarField::farFieldTemperatureDSTFvPatchScalarField
(
    const farFieldTemperatureDSTFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    UInf_(ptf.UInf_),
    pInf_(ptf.pInf_),
    TInf_(ptf.TInf_)
{}


farFieldTemperatureDSTFvPatchScalarField::farFieldTemperatureDSTFvPatchScalarField
(
    const farFieldTemperatureDSTFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    UInf_(ptf.UInf_),
    pInf_(ptf.pInf_),
    TInf_(ptf.TInf_)
{}


farFieldTemperatureDSTFvPatchScalarField::farFieldTemperatureDSTFvPatchScalarField
(
    const farFieldTemperatureDSTFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    UInf_(ptf.UInf_),
    pInf_(ptf.pInf_),
    TInf_(ptf.TInf_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void farFieldTemperatureDSTFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    label patchIndex = patch().index();
    const basicThermo& thermo = db().lookupObject<basicThermo>("thermophysicalProperties");
    scalarField phi = db().lookupObject<surfaceScalarField>("phi").boundaryField()[patchIndex];
    scalarField rhoP = db().lookupObject<volScalarField>("rho").boundaryField()[patchIndex];
    scalarField rhoI = db().lookupObject<volScalarField>("rho").boundaryField()[patchIndex].patchInternalField();
    vectorField UI = db().lookupObject<volVectorField>("U").boundaryField()[patchIndex].patchInternalField();
    scalarField pP = db().lookupObject<volScalarField>("p").boundaryField()[patchIndex];
    scalarField pI = db().lookupObject<volScalarField>("p").boundaryField()[patchIndex].patchInternalField();
    
    scalarField c = 
        sqrt
        (
            thermo.Cp()().boundaryField()[patchIndex]
          / thermo.Cv()().boundaryField()[patchIndex]
          / thermo.psi().boundaryField()[patchIndex]
        );
    scalarField cI = 
        sqrt
        (
            thermo.Cp()().boundaryField()[patchIndex].patchInternalField()
          / thermo.Cv()().boundaryField()[patchIndex].patchInternalField()
          / thermo.psi().boundaryField()[patchIndex].patchInternalField()
        );
    
    scalarField mach = phi/rhoP/patch().magSf()/c;


    scalarField supersonicInflow(mach.size(),0);
    scalarField supersonicOutflow(mach.size(),0);
        
    DST::neg0(supersonicInflow,(mach+1)());
    DST::pos0(supersonicOutflow,(mach-1)());
    

    scalarField subsonicInflow(mach.size(),0);
    scalarField subsonicOutflow(mach.size(),0);
    
    DST::neg(subsonicInflow, mach);
    DST::pos0(subsonicOutflow, mach);
    
    subsonicInflow -= supersonicInflow;
    subsonicOutflow -= supersonicOutflow;
    
    
    vectorField n(patch().nf());
    scalarField rho0c0 = rhoI*cI;
    scalarField pb = 0.5*
        (
            pInf_
          + pI
          - rho0c0
          * (
                n & (UInf_ - UI)
            )
        );
    
    scalarField RI = pI/rhoI/patchInternalField();
    scalarField rhoa = pInf_/TInf_/RI;

    scalarField subsInValue = pb/(rhoa+(pb-pInf_)/sqr(cI))/RI;
    scalarField subsOutValueFraction = pb/(pI+RI*patchInternalField()*(pb-pI)/sqr(cI));

    this->refValue() = 
        (
            supersonicInflow * TInf_
          + subsonicInflow * subsInValue
        );
    
    this->valueFraction() = 1 - supersonicOutflow - subsOutValueFraction*subsonicOutflow;

    mixedFvPatchScalarField::updateCoeffs();
}


void farFieldTemperatureDSTFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    os.writeKeyword("UInf") << UInf_ << token::END_STATEMENT << nl;
    os.writeKeyword("pInf") << pInf_ << token::END_STATEMENT << nl;
    os.writeKeyword("TInf") << TInf_ << token::END_STATEMENT << nl;

    this->writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, farFieldTemperatureDSTFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
