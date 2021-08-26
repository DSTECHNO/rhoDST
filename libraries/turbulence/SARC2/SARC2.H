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
    Foam::SARC2

Description
    Simpler Rotation-Curvature Correction method is applied to Spalart-Allmaras
    turbulence model

SourceFiles
    SARC2.C

References
-----------
    Z. Qiang, Y. Yong, A New, Simpler Rotation/Curvature Correction Method for the Spalart-Allmaras 
    Turbulence Model, Chinese Journal of Aeronautics Vol. 26, 2, 2013, pp. 326-333.   
 
\*---------------------------------------------------------------------------*/

#ifndef compressibleSARC2_H
#define compressibleSARC2_H

#include "RASModel.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                      Class SARC2 Declaration
\*---------------------------------------------------------------------------*/

class SARC2
:
    public RASModel
{
    // Private data

        // Model coefficients

            dimensionedScalar sigmaNut_;
            dimensionedScalar kappa_;
            dimensionedScalar Prt_;

            dimensionedScalar Cb1_;
            dimensionedScalar Cb2_;
            dimensionedScalar Cw1_;
            dimensionedScalar Cw2_;
            dimensionedScalar Cw3_;
            dimensionedScalar Cv1_;
            dimensionedScalar Cv2_;

        dimensionedScalar Cr1_;
            dimensionedScalar Cr2_;
            dimensionedScalar Cr3_;

        Switch RCMCorrection_;


        // Fields

            volScalarField nuTilda_;
            volScalarField mut_;
            volScalarField alphat_;


        //- Wall distance
        wallDist d_;


    // Private member functions

        tmp<volScalarField> chi() const;
        tmp<volScalarField> fv1(const volScalarField& chi) const;
        tmp<volScalarField> fv2
        (
            const volScalarField& chi,
            const volScalarField& fv1
        ) const;
        tmp<volScalarField> fv3
        (
            const volScalarField& chi,
            const volScalarField& fv1
        ) const;
        tmp<volScalarField> fw(const volScalarField& Stilda) const;

    tmp<volScalarField> fr1(const volTensorField& gradU) const;

public:

    //- Runtime type information
    TypeName("SARC2");


    // Constructors

        //- Construct from components
        SARC2
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const basicThermo& thermophysicalModel,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~SARC2()
    {}


    // Member Functions

        //- Return the effective diffusivity for nuTilda
        tmp<volScalarField> DnuTildaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DnuTildaEff",
                    rho_*nuTilda_/sigmaNut_ + mu()
                )
            );
        }

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> mut() const
        {
            return mut_;
        }

        //- Return the turbulence thermal diffusivity
        virtual tmp<volScalarField> alphat() const
        {
            return alphat_;
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "k",
                        runTime_.timeName(),
                        mesh_
                    ),
                    mesh_,
                    dimensionedScalar("0", dimensionSet(0, 2, -2, 0, 0), 0)
                )
            );
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "epslion",
                        runTime_.timeName(),
                        mesh_
                    ),
                    mesh_,
                    dimensionedScalar("0", dimensionSet(0, 2, -3, 0, 0), 0)
                )
            );
        }

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devRhoReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff() const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
