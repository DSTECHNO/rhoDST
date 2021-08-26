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

\*-------------------------------------------------------------------------------------*/

#include "steadyStateDST.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
steadyStateDST<Type>::fvcDdt
(
    const dimensioned<Type>& dt
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+dt.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    return tmp<GeometricField<Type, fvPatchField, volMesh> >
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            ddtIOobject,
            mesh(),
            dimensioned<Type>
            (
                "0",
                dt.dimensions()/dimTime,
                pTraits<Type>::zero
            ),
            calculatedFvPatchField<Type>::typeName
        )
    );
}

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
steadyStateDST<Type>::fvcDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const volScalarField& rDeltaT = mesh().template lookupObject<volScalarField>("rDeltaTDST");

    IOobject ddtIOobject
    (
        "ddt("+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    return tmp<GeometricField<Type, fvPatchField, volMesh> >
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            ddtIOobject,
            rDeltaT*(vf - vf.oldTime())
        )
    );
}

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
steadyStateDST<Type>::fvcDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const volScalarField& rDeltaT = mesh().template lookupObject<volScalarField>("rDeltaTDST");

    IOobject ddtIOobject
    (
        "ddt("+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    return tmp<GeometricField<Type, fvPatchField, volMesh> >
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            ddtIOobject,
            rDeltaT*rho*(vf - vf.oldTime())
        )
    );
}

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
steadyStateDST<Type>::fvcDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const volScalarField& rDeltaT = mesh().template lookupObject<volScalarField>("rDeltaTDST");

    IOobject ddtIOobject
    (
        "ddt("+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    return tmp<GeometricField<Type, fvPatchField, volMesh> >
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            ddtIOobject,
            rDeltaT*(rho*vf - rho.oldTime()*vf.oldTime())
        )
    );
}

template<class Type>
tmp<fvMatrix<Type> >
steadyStateDST<Type>::fvmDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*dimVol/dimTime
        )
    );

    fvMatrix<Type>& fvm = tfvm();

    const volScalarField& rDeltaT = mesh().template lookupObject<volScalarField>("rDeltaTDST");

    fvm.diag() = rDeltaT.internalField()*mesh().V();

    fvm.source() = rDeltaT*vf.oldTime().internalField()*mesh().V();

    return tfvm;
}

template<class Type>
tmp<fvMatrix<Type> >
steadyStateDST<Type>::fvmDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm();

    const volScalarField& rDeltaT = mesh().template lookupObject<volScalarField>("rDeltaTDST");

    fvm.diag() = rDeltaT.internalField()*rho.value()*mesh().V();

    fvm.source() = rDeltaT.internalField()
        *rho.value()*vf.oldTime().internalField()*mesh().V();

    return tfvm;
}

template<class Type>
tmp<fvMatrix<Type> >
steadyStateDST<Type>::fvmDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm();

    const volScalarField& rDeltaT = mesh().template lookupObject<volScalarField>("rDeltaTDST");

    fvm.diag() = rDeltaT.internalField()*rho.internalField()*mesh().V();

    fvm.source() = rDeltaT.internalField()
        *rho.oldTime().internalField()
        *vf.oldTime().internalField()*mesh().V();

    return tfvm;
}

template<class Type>
tmp<typename steadyStateDST<Type>::fluxFieldType>
steadyStateDST<Type>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phiAbs
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddtPhiCorr(" + rA.name() + ',' + U.name() + ',' + phiAbs.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    tmp<fluxFieldType> phiCorr =
        phiAbs.oldTime() - (fvc::interpolate(U.oldTime()) & mesh().Sf());

    return tmp<fluxFieldType>
    (
        new fluxFieldType
        (
            ddtIOobject,
            this->fvcDdtPhiCoeff(U.oldTime(), phiAbs.oldTime(), phiCorr())
           *fvc::interpolate(rDeltaT*rA)*phiCorr
        )
    );
}

template<class Type>
tmp<typename steadyStateDST<Type>::fluxFieldType>
steadyStateDST<Type>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phiAbs
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddtPhiCorr("
      + rA.name() + ','
      + rho.name() + ','
      + U.name() + ','
      + phiAbs.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    if
    (
        U.dimensions() == dimVelocity
     && phiAbs.dimensions() == dimVelocity*dimArea
    )
    {
        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                rDeltaT
               *this->fvcDdtPhiCoeff(U.oldTime(), phiAbs.oldTime())
               *(
                   fvc::interpolate(rA*rho.oldTime())*phiAbs.oldTime()
                 - (fvc::interpolate(rA*rho.oldTime()*U.oldTime())
                  & mesh().Sf())
                )
            )
        );
    }
    else if
    (
        U.dimensions() == dimVelocity
     && phiAbs.dimensions() == rho.dimensions()*dimVelocity*dimArea
    )
    {
        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                rDeltaT
               *this->fvcDdtPhiCoeff
                (
                    U.oldTime(),
                    phiAbs.oldTime()/fvc::interpolate(rho.oldTime())
                )
               *(
                   fvc::interpolate(rA*rho.oldTime())
                  *phiAbs.oldTime()/fvc::interpolate(rho.oldTime())
                 - (
                       fvc::interpolate
                       (
                           rA*rho.oldTime()*U.oldTime()
                       ) & mesh().Sf()
                   )
                )
            )
        );
    }
    else if
    (
        U.dimensions() == rho.dimensions()*dimVelocity
     && phiAbs.dimensions() == rho.dimensions()*dimVelocity*dimArea
    )
    {
        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                rDeltaT
               *this->fvcDdtPhiCoeff
                (
                    rho.oldTime(),
                    U.oldTime(),
                    phiAbs.oldTime()
                )
               *(
                   fvc::interpolate(rA)*phiAbs.oldTime()
                 - (fvc::interpolate(rA*U.oldTime()) & mesh().Sf())
                )
            )
        );
    }
    else
    {
        FatalErrorIn
        (
            "steadyStateDST<Type>::fvcDdtPhiCorr"
        )   << "dimensions of phiAbs are not correct"
            << abort(FatalError);

        return fluxFieldType::null();
    }
}

template<class Type>
tmp<typename steadyStateDST<Type>::fluxFieldType>
steadyStateDST<Type>::fvcDdtConsistentPhiCorr
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& faceU,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const surfaceScalarField& rAUf
)
{
    tmp<fluxFieldType> toldTimeFlux =
        (mesh().Sf() & faceU.oldTime())*rAUf/mesh().time().deltaT();

    if (mesh().moving())
    {
        // Mesh is moving, need to take into account the ratio between old and
        // current cell volumes
        volScalarField V0ByV
        (
            IOobject
            (
                "V0ByV",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("one", dimless, 1.0),
            zeroGradientFvPatchScalarField::typeName
        );
        V0ByV.internalField() = mesh().V0()/mesh().V();
        V0ByV.correctBoundaryConditions();

        // Correct the flux with interpolated volume ratio
        toldTimeFlux() *= fvc::interpolate(V0ByV);
    }

    return toldTimeFlux;
}

template<class Type>
tmp<surfaceScalarField> steadyStateDST<Type>::meshPhi
(
    const GeometricField<Type, fvPatchField, volMesh>&
)
{
    return mesh().phi();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
