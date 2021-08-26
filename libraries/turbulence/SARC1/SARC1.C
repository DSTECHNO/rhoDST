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

\*---------------------------------------------------------------------------*/

#include "SARC1.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcDdt.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(SARC1, 0);
addToRunTimeSelectionTable(RASModel, SARC1, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> SARC1::chi() const
{
    return rho_*nuTilda_/mu();
}


tmp<volScalarField> SARC1::fv1(const volScalarField& chi) const
{
    volScalarField chi3 = pow3(chi);
    return chi3/(chi3 + pow3(Cv1_));
}


tmp<volScalarField> SARC1::fv2
(
    const volScalarField& chi,
    const volScalarField& fv1
) const
{
    return 1.0/pow3(scalar(1) + chi/Cv2_);
}


tmp<volScalarField> SARC1::fv3
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


tmp<volScalarField> SARC1::fw(const volScalarField& Stilda) const
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
                      SARC1 modification
\*---------------------------------------------------------------------------*/
tmp<volScalarField> SARC1::fr1(const volTensorField& S, const volTensorField& W) const
{
    
    if (spalartShurCorrection_)
    {
        volScalarField sqrS(2.0*magSqr(S));
        volScalarField sqrW(2.0*magSqr(W));
        volScalarField sqrD(0.5*(sqrS + sqrW));

        volScalarField rStar(sqrt(sqrS/max(sqrW, dimensionedScalar("SMALL", sqrW.dimensions(), SMALL))));

        volScalarField Wxx(W.component(tensor::XX));
        volScalarField Wxy(W.component(tensor::XY));
        volScalarField Wxz(W.component(tensor::XZ));
        volScalarField Wyx(W.component(tensor::YX));
        volScalarField Wyy(W.component(tensor::YY));
        volScalarField Wyz(W.component(tensor::YZ));
        volScalarField Wzx(W.component(tensor::ZX));
        volScalarField Wzy(W.component(tensor::ZY));
        volScalarField Wzz(W.component(tensor::ZZ));

    
        
        volScalarField rTilda
        (
            2.0/(rho_*sqr(max(sqrD, dimensionedScalar("SMALL", sqrD.dimensions(), SMALL))))
               *(
                   (Wxx*Sxx_ + Wxy*Sxy_ + Wxz*Sxz_)
                      *(fvc::ddt(rho_,Sxx_)+fvc::div(phi_, Sxx_)) //i,j=1,1
                 + (Wxx*Syx_ + Wxy*Syy_ + Wxz*Syz_)
                      *( fvc::ddt(rho_,Sxy_)+fvc::div(phi_, Sxy_)) //i,j=1,2
                 + (Wxx*Szx_ + Wxy*Szy_ + Wxz*Szz_)
                      *(fvc::ddt(rho_,Sxz_)+fvc::div(phi_, Sxz_)) //i,j=1,3
                 + (Wyx*Sxx_ + Wyy*Sxy_ + Wyz*Sxz_)
                      *(fvc::ddt(rho_,Syx_)+fvc::div(phi_, Syx_)) //i,j=2,1
                 + (Wyx*Syx_ + Wyy*Syy_ + Wyz*Syz_)
                      *(fvc::ddt(rho_,Syy_)+fvc::div(phi_, Syy_)) //i,j=2,2
                 + (Wyx*Szx_ + Wyy*Szy_ + Wyz*Szz_)
                      *( fvc::ddt(rho_,Syz_)+fvc::div(phi_, Syz_)) //i,j=2,3
                 + (Wzx*Sxx_ + Wzy*Sxy_ + Wzz*Sxz_)
                      *(fvc::ddt(rho_,Szx_)+fvc::div(phi_, Szx_)) //i,j=3,1
                 + (Wzx*Syx_ + Wzy*Syy_ + Wzz*Syz_)
                      *( fvc::ddt(rho_,Szy_)+fvc::div(phi_, Szy_)) //i,j=3,2
                 + (Wzx*Szx_ + Wzy*Szy_ + Wzz*Szz_)
                      *(fvc::ddt(rho_,Szz_)+fvc::div(phi_, Szz_)) //i,j=3,3
                )
        );

        return 
            (1 + Cr1_)*2.0*rStar/(1.0 + rStar)*(1.0 - Cr3_*atan(Cr2_*rTilda))
             - Cr1_;
    }
    else
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "fr1",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("fr1", dimless, 1),
                zeroGradientFvPatchScalarField::typeName
            )
        );
    }
}
/*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

SARC1::SARC1
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
                      SARC1 modification
\*---------------------------------------------------------------------------*/
    Cr1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cr1",
            coeffDict_,
            1.0
        )
    ),
    Cr2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cr2",
            coeffDict_,
            12.0
        )
    ),
    Cr3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cr3",
            coeffDict_,
            1.0
        )
    ),

    spalartShurCorrection_
    (
        coeffDict_.lookupOrDefault("spalartShurCorrection", true)
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
/*---------------------------------------------------------------------------*\
                      SARC1 modification
\*---------------------------------------------------------------------------*/
    Sxx_
    (
        IOobject
        (
            "Sxx",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Sxx", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    Sxy_
    (
        IOobject
        (
            "Sxy",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Sxy", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    Sxz_
    (
        IOobject
        (
            "Sxz",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Sxz", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    Syx_
    (
        IOobject
        (
            "Syx",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Syx", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    Syy_
    (
        IOobject
        (
            "Syy",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Syy", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    Syz_
    (
        IOobject
        (
            "Syz",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Syz", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    Szx_
    (
        IOobject
        (
            "Szx",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Szx", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    Szy_
    (
        IOobject
        (
            "Szy",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Szy", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    Szz_
    (
        IOobject
        (
            "Szz",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Szz", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0),
        zeroGradientFvPatchScalarField::typeName
    ),

/*---------------------------------------------------------------------------*/

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

    if (spalartShurCorrection_)
    {
        Info<< "    Employing Spalart-Shur Rotation/Curvature correction model for the turbulence" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> SARC1::R() const
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


tmp<volSymmTensorField> SARC1::devRhoReff() const
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


tmp<fvVectorMatrix> SARC1::divDevRhoReff() const
{
    volScalarField muEff_ = muEff();

    return
    (
      - fvm::laplacian(muEff_, U_)
      - fvc::div(muEff_*dev2(T(fvc::grad(U_))))
    );
}


bool SARC1::read()
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
                      SARC1 modification
\*---------------------------------------------------------------------------*/
    Cr1_.readIfPresent(coeffDict());
        Cr2_.readIfPresent(coeffDict());
        Cr3_.readIfPresent(coeffDict());

    spalartShurCorrection_.readIfPresent
        (
            "spalartShurCorrection", coeffDict()
        );
/*---------------------------------------------------------------------------*/

        return true;
    }
    else
    {
        return false;
    }
}


void SARC1::correct()
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

    volScalarField Stilda =
        fv3(chi, fv1)*::sqrt(2.0)*mag(skew(fvc::grad(U_)))
      + fv2(chi, fv1)*nuTilda_/sqr(kappa_*d_);
/*---------------------------------------------------------------------------*\
                      SARC1 modification
\*---------------------------------------------------------------------------*/
    volTensorField gradU(fvc::grad(U_));
    volTensorField S((gradU+gradU.T())/2);
    volTensorField W(skew(gradU));

    Sxx_ = S.component(tensor::XX);
    Sxy_ = S.component(tensor::XY);
    Sxz_ = S.component(tensor::XZ);
    Syx_ = S.component(tensor::YX);
    Syy_ = S.component(tensor::YY);
    Syz_ = S.component(tensor::YZ);
    Szx_ = S.component(tensor::ZX);
    Szy_ = S.component(tensor::ZY);
    Szz_ = S.component(tensor::ZZ);
/*---------------------------------------------------------------------------*/
    tmp<fvScalarMatrix> nuTildaEqn
    (
        fvm::ddt(rho_, nuTilda_)
      + fvm::div(phi_, nuTilda_)
      - fvm::laplacian(DnuTildaEff(), nuTilda_)
      - Cb2_/sigmaNut_*rho_*magSqr(fvc::grad(nuTilda_))
     ==
        Cb1_*rho_*Stilda*nuTilda_*fr1(S, W)
      - fvm::Sp(Cw1_*fw(Stilda)*nuTilda_*rho_/sqr(d_), nuTilda_)
    );
    nuTildaEqn().relax();
    solve(nuTildaEqn);
    bound(nuTilda_, dimensionedScalar("0", nuTilda_.dimensions(), 0.0));
    nuTilda_.correctBoundaryConditions();

    // Re-calculate viscosity
    mut_.internalField() = fv1*nuTilda_.internalField()*rho_.internalField();
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
