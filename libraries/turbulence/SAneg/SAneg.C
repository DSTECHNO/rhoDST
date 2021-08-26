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

#include "SAneg.H"
#include "addToRunTimeSelectionTable.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(SAneg, 0);
addToRunTimeSelectionTable(RASModel, SAneg, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> SAneg::chi() const
{
    return rho_*nuTilda_/mu();
}


tmp<volScalarField> SAneg::fv1(const volScalarField& chi) const
{
    volScalarField chi3 = pow3(chi);
    return chi3/(chi3 + pow3(Cv1_));
}


tmp<volScalarField> SAneg::fv2
(
    const volScalarField& chi,
    const volScalarField& fv1
) const
{
    return 1.0/pow3(scalar(1) + chi/Cv2_);
}


tmp<volScalarField> SAneg::fv3
(
    const volScalarField& chi,
    const volScalarField& fv1
) const
{
    volScalarField chiByCv2 = (1/Cv2_)*chi;

    return
        (scalar(1) + chi*fv1)
       *(1/Cv2_)
       *(3*(scalar(1) + chiByCv2) + sqr(chiByCv2))
       /pow3(scalar(1) + chiByCv2);
}


tmp<volScalarField> SAneg::fw(const volScalarField& Stilda) const
{
    volScalarField r = min
    (
        nuTilda_
       /(
           max(Stilda, dimensionedScalar("SMALL", Stilda.dimensions(), SMALL))
           *sqr(kappa_*d_)
        ),
        scalar(10.0)
    );
    r.boundaryField() == 0.0;

    volScalarField g = r + Cw2_*(pow6(r) - r);

    return g*pow((1.0 + pow6(Cw3_))/(pow6(g) + pow6(Cw3_)), 1.0/6.0);
}

/*---------------------------------------------------------------------------*\
                      SAneg modification
\*---------------------------------------------------------------------------*/
tmp<volScalarField> SAneg::ft2(const volScalarField& chi) const
{
    if(negativeNuTilda_)
    {
        return Ct3_*exp(-Ct4_*sqr(chi));
    }
    else
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "ft2",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("ft2", dimless, 0),
                zeroGradientFvPatchScalarField::typeName
            )
        );
    }
}
/*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

SAneg::SAneg
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const basicThermo& thermophysicalModel,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, rho, U, phi, thermophysicalModel, turbulenceModelName),

    sigmaNut_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaNut",
            coeffDict_,
            0.66666
        )
    ),
    kappa_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "kappa",
            coeffDict_,
            0.41
        )
    ),
    Prt_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Prt",
            coeffDict_,
            1.0
        )
    ),

    Cb1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Cb1",
            coeffDict_,
            0.1355
        )
    ),
    Cb2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Cb2",
            coeffDict_,
            0.622
        )
    ),
    Cw1_(Cb1_/sqr(kappa_) + (1.0 + Cb2_)/sigmaNut_),
    Cw2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Cw2",
            coeffDict_,
            0.3
        )
    ),
    Cw3_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Cw3",
            coeffDict_,
            2.0
        )
    ),
    Cv1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Cv1",
            coeffDict_,
            7.1
        )
    ),
    Cv2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Cv2",
            coeffDict_,
            5.0
        )
    ),
/*---------------------------------------------------------------------------*\
                      SAneg modification
\*---------------------------------------------------------------------------*/
     Cn1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cn1",
            coeffDict_,
            16.0
        )
    ),
    Cn2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cn2",
            coeffDict_,
            0.7
        )
    ),
    Cn3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cn3",
            coeffDict_,
            0.9
        )
    ),
    Ct3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ct3",
            coeffDict_,
            1.2
        )
    ),
    Ct4_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ct4",
            coeffDict_,
            0.5
        )
    ),
    StildaModification_
    (
        coeffDict_.lookupOrDefault("StildaModification", true)
    ),

    negativeNuTilda_
    (
        coeffDict_.lookupOrDefault("negativeNuTilda", true)
    ),
/*---------------------------------------------------------------------------*/

    nuTilda_
    (
        IOobject
        (
            "nuTilda",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    mut_
    (
        IOobject
        (
            "mut",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    alphat_
    (
        IOobject
        (
            "alphat",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateAlphat("alphat", mesh_)
    ),

    d_(mesh_)
{
    alphat_ = mut_/Prt_;
    alphat_.correctBoundaryConditions();

    printCoeffs();

/*---------------------------------------------------------------------------*\
                      SAneg modification
\*---------------------------------------------------------------------------*/
    if (StildaModification_)
    {
        Info<< "    Enabling new Stilda modification" << endl;
    }
    if (negativeNuTilda_)
    {
        Info<< "    Enabling negative nuTilda" << endl;
    }
}
/*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                      SAneg modification
\*---------------------------------------------------------------------------*/
tmp<volScalarField> SAneg::DnuTildaEff(const volScalarField& chi) const
{
    volScalarField pow3chi(pow(chi, 3));
    volScalarField fn
    (
        pos(chi) + neg(chi)*(Cn1_ + pow3chi)/(Cn1_ - pow3chi)
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "DnuTildaEff",
            (
                (rho_*nuTilda_*fn + mu())/sigmaNut_
            )
        )
    );
}
/*---------------------------------------------------------------------------*/

tmp<volSymmTensorField> SAneg::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k() - (mut_/rho_)*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<volSymmTensorField> SAneg::devRhoReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -muEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> SAneg::divDevRhoReff() const
{
    volScalarField muEff_ = muEff();

    return
    (
      - fvm::laplacian(muEff_, U_)
      - fvc::div(muEff_*dev2(T(fvc::grad(U_))))
    );
}


bool SAneg::read()
{
    if (RASModel::read())
    {
        sigmaNut_.readIfPresent(coeffDict());
        kappa_.readIfPresent(coeffDict());
        Prt_.readIfPresent(coeffDict());

        Cb1_.readIfPresent(coeffDict());
        Cb2_.readIfPresent(coeffDict());
        Cw1_ = Cb1_/sqr(kappa_) + (1.0 + Cb2_)/sigmaNut_;
        Cw2_.readIfPresent(coeffDict());
        Cw3_.readIfPresent(coeffDict());
        Cv1_.readIfPresent(coeffDict());
        Cv2_.readIfPresent(coeffDict());
/*---------------------------------------------------------------------------*\
                      SAneg modification
\*---------------------------------------------------------------------------*/
        Cn1_.readIfPresent(coeffDict());
        Cn2_.readIfPresent(coeffDict());
        Cn3_.readIfPresent(coeffDict());
        Ct3_.readIfPresent(coeffDict());
        Ct4_.readIfPresent(coeffDict());

        StildaModification_.readIfPresent
        (
            "StildaModification", coeffDict()
        );
        negativeNuTilda_.readIfPresent
        (
            "negativeNuTilda", coeffDict()
        );
/*---------------------------------------------------------------------------*/
        return true;
    }
    else
    {
        return false;
    }
}


void SAneg::correct()
{
    // Bound in case of topological change
    // HJ, 22/Aug/2007
    if (mesh_.changing())
    {
        bound(nuTilda_, dimensionedScalar("0", nuTilda_.dimensions(), 0.0));
    }

    if (!turbulence_)
    {
        // Re-calculate viscosity
        mut_ = rho_*nuTilda_*fv1(chi());
        mut_.correctBoundaryConditions();

        // Re-calculate thermal diffusivity
        alphat_ = mut_/Prt_;
        alphat_.correctBoundaryConditions();

        return;
    }

    RASModel::correct();

    if (mesh_.changing())
    {
        d_.correct();
    }

    volScalarField chi = this->chi();
    volScalarField fv1 = this->fv1(chi);
/*---------------------------------------------------------------------------*\
                      SAneg modification
\*---------------------------------------------------------------------------*/
    volTensorField gradU(fvc::grad(U_));
    volSymmTensorField S(symm(gradU));
    volTensorField W(skew(gradU));
    volScalarField Omega(sqrt(2.0)*mag(W));
    
    volScalarField Stilda(
        IOobject
        (
            "Stilda",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Stilda", Omega.dimensions(), 0.0)
    );
    if (StildaModification_)
    {
        volScalarField Sbar
        (
            fv2(chi, fv1)*nuTilda_/sqr(kappa_*d_)
        );
        Stilda = Omega
               + pos(Cn2_*Omega + Sbar)*Sbar
               + neg(Cn2_*Omega + Sbar)
               *(Omega*(sqr(Cn2_)*Omega + Cn3_*Sbar))
               /max(((Cn3_ - 2.0*Cn2_)*Omega - Sbar),
                    dimensionedScalar("SMALL", Omega.dimensions(), SMALL));
    }
    else
    {
        Stilda = fv3(chi, fv1)*::sqrt(2.0)*mag(W)
               + fv2(chi, fv1)*nuTilda_/sqr(kappa_*d_);
    }
/*---------------------------------------------------------------------------*/
    
    tmp<fvScalarMatrix> nuTildaEqn
    (
        fvm::ddt(rho_, nuTilda_)
      + fvm::div(phi_, nuTilda_)
/*---------------------------------------------------------------------------*\
                      SAneg modification
\*---------------------------------------------------------------------------*/
      - fvm::laplacian(DnuTildaEff(chi), nuTilda_)
/*---------------------------------------------------------------------------*/
      - Cb2_/sigmaNut_*rho_*magSqr(fvc::grad(nuTilda_))
     ==
       pos(nuTilda_)
       *(
        Cb1_*(1.0 - ft2(chi))*rho_*Stilda*nuTilda_
              - fvm::Sp((Cw1_*fw(Stilda) -
            Cb1_/sqr(kappa_)*ft2(chi))*nuTilda_*rho_/sqr(d_), nuTilda_)
        )
       + neg(nuTilda_)
       *(
            Cb1_*(1.0 - Ct3_)*rho_*Omega*nuTilda_
          + fvm::Sp(Cw1_*nuTilda_*rho_/sqr(d_), nuTilda_)
        )
    );

    nuTildaEqn().relax();
    solve(nuTildaEqn);
    if(!negativeNuTilda_)
    {
        // bound nuTilda when negativeNuTilda = false
        bound(nuTilda_, dimensionedScalar("0", nuTilda_.dimensions(), 0.0));
    }
    nuTilda_.correctBoundaryConditions();

    // Re-calculate viscosity
    mut_.internalField() = mag(fv1)*nuTilda_.internalField()*rho_.internalField();
    if(negativeNuTilda_)
    {
        bound(mut_, dimensionedScalar("0", mut_.dimensions(), 0.0));
    }
    mut_.correctBoundaryConditions();

    // Re-calculate thermal diffusivity
    alphat_ = mut_/Prt_;
    alphat_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
