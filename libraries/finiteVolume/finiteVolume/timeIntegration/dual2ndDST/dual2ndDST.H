/*-----------------------------------------------------------------------------------------*\
| =========                 |                                                               |
| \\      /  F ield         | rhoDST: Block-coupled solver for compressible flows           |
|  \\    /   O peration     |								    |
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
    rhoDST::dual2ndDST

Description
    Second-order dual time ddt scheme. Global and local deltaT is used for 
    outer and inner iterations respectively.

SourceFiles
    dual2ndDST.C

\*---------------------------------------------------------------------------------------*/

#ifndef dual2ndDST_H
#define dual2ndDST_H

#include "ddtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class dual2ndDST Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class dual2ndDST
:
    public fv::ddtScheme<Type>
{
    // Private Member Functions

        //- Return the current time-step
            scalar deltaT_() const;

        //- Return the previous time-step
            scalar deltaT0_() const;

        //- Return the previous time-step or GREAT if the old timestep field
        //  wasn't available in which case Euler ddt is used
            template<class GeoField>
            scalar deltaT0_(const GeoField&) const;

        //- Disallow default bitwise copy construct
            dual2ndDST(const dual2ndDST&);

        //- Disallow default bitwise assignment
            void operator=(const dual2ndDST&);

public:

    //- Runtime type information
        TypeName("dual2ndDST");


    // Constructors

        //- Construct from mesh
            dual2ndDST(const fvMesh& mesh)
            :
                ddtScheme<Type>(mesh)
            {}

        //- Construct from mesh and Istream
            dual2ndDST(const fvMesh& mesh, Istream& is)
            :
                ddtScheme<Type>(mesh, is)
            {}


    // Member Functions

        //- Return mesh reference
            const fvMesh& mesh() const
            {
                return fv::ddtScheme<Type>::mesh();
            }

            tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt
            (
                const dimensioned<Type>&
            );

            tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt
            (
                const GeometricField<Type, fvPatchField, volMesh>&
            );

            tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt
            (
                const dimensionedScalar&,
                const GeometricField<Type, fvPatchField, volMesh>&
            );

            tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt
            (
                const volScalarField&,
                const GeometricField<Type, fvPatchField, volMesh>&
            );

            tmp<fvMatrix<Type> > fvmDdt
            (
                const GeometricField<Type, fvPatchField, volMesh>&
            );

            tmp<fvMatrix<Type> > fvmDdt
            (
                const dimensionedScalar&,
                const GeometricField<Type, fvPatchField, volMesh>&
            );

            tmp<fvMatrix<Type> > fvmDdt
            (
                const volScalarField&,
                const GeometricField<Type, fvPatchField, volMesh>&
            );

            typedef typename ddtScheme<Type>::fluxFieldType fluxFieldType;

            tmp<fluxFieldType> fvcDdtPhiCorr
            (
                const volScalarField& rA,
                const GeometricField<Type, fvPatchField, volMesh>& U,
                const fluxFieldType& phi
            );

            tmp<fluxFieldType> fvcDdtPhiCorr
            (
                const volScalarField& rA,
                const volScalarField& rho,
                const GeometricField<Type, fvPatchField, volMesh>& U,
                const fluxFieldType& phi
            );


            // Member functions for the new time consistent formulation
            tmp<fluxFieldType> fvcDdtConsistentPhiCorr
            (
                const GeometricField<Type, fvsPatchField, surfaceMesh>& faceU,
                const GeometricField<Type, fvPatchField, volMesh>& U,
                const surfaceScalarField& rAUf
            );


            tmp<surfaceScalarField> meshPhi
            (
                const GeometricField<Type, fvPatchField, volMesh>&
            );
};


template<>
tmp<surfaceScalarField> dual2ndDST<scalar>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const volScalarField& U,
    const surfaceScalarField& phi
);


template<>
tmp<surfaceScalarField> dual2ndDST<scalar>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const volScalarField& rho,
    const volScalarField& U,
    const surfaceScalarField& phi
);


template<>
tmp<surfaceScalarField> dual2ndDST<scalar>::fvcDdtConsistentPhiCorr
(
    const surfaceScalarField& faceU,
    const volScalarField& U,
    const surfaceScalarField& rAUf
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "dual2ndDST.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
