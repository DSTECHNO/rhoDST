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

#include "DST.H"
#include "timeDST.H"


namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

timeDST::timeDST
(
    const argList& args
)
:
    Time(Foam::Time::controlDictName, args)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

timeDST::~timeDST()
{}

// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

void timeDST::initialize
(
    basicThermo& thermo,
    const IOdictDST& dict,
    compressible::turbulenceModel& turbulence,
    volVectorField& U
)
{
    thermoPtr = &thermo;
    turbulencePtr = &turbulence;
    UPtr = &U;
    const fvMesh& mesh = UPtr->mesh();

    dictPtr = &dict;

    steadyState = false;

    word ddtSc(mesh.schemesDict().ddtScheme("default"));

    if(ddtSc == "steadyStateDST" || ddtSc == "steadyState")
        steadyState = true;

    if(steadyState)
    {
        rhoOldPtr = autoPtr<volScalarField>
        (   new volScalarField
            (
                IOobject
                (
                    "rhoOld",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                thermoPtr->rho()
            )
        );
        
        eOldPtr = autoPtr<volScalarField> 
        (   new volScalarField
            (
                IOobject
                (
                    "eOld",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                thermoPtr->h()*thermoPtr->Cv()/thermoPtr->Cp()
            )
        );
    }

    neighbouring = dictPtr->neighbouring();

    if(!neighbouring)
    {
        // courant is defined at each cell center as described by Blazek(2001)
        const surfaceVectorField& Sf = mesh.Sf();
        vectorField magSf = cmptMag(mesh.Sf().internalField());

        dS = autoPtr<volVectorField>
        (
            new volVectorField
            (
                IOobject
                (
                    "dS",
                    UPtr->time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedVector("dS",vector(0,0,0))
            )
        );
        
        const labelList& owner = mesh.owner();
        const labelList& neighbour=  mesh.neighbour();
        forAll(magSf, faceI)
        {
            dS()[owner[faceI]] += magSf[faceI];
            dS()[neighbour[faceI]] += magSf[faceI];
        }
        forAll(Sf.boundaryField(), patchI)
        {
            const labelUList& fc = mesh.boundary()[patchI].faceCells();
            vectorField Sfp = cmptMag(Sf.boundaryField()[patchI]);
            forAll(fc, cellI)
            {
                dS()[fc[cellI]] += Sfp[cellI];
            }
        }
        dS() *= 0.5;
    }

    residualRatio = 1;

    maxCo = dictPtr->maxCo();
    minCo = dictPtr->minCo();

    courant = dictPtr->initialCo();

    courantPtr = autoPtr<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "courant",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("courant",dimless,courant)
        )
    );

    innerIter = 0;
    maxInnerIter = dictPtr->lookupOrDefault<int>("maxIter", 10);

    rDeltaTDSTPtr = autoPtr<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "rDeltaTDST",
                U.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("rDeltaTDST",dimless/dimTime,0)
        )
    );
    
    if(steadyState)
    {
        setSteadyStateDeltaT();
    }
    else
    {
        setTransientDeltaT();
        courantPtr->internalField() = courant;
        setSteadyStateDeltaT();
    }
}

void timeDST::setSteadyStateDeltaT()
{
    const fvMesh& mesh = UPtr->mesh();
    if(!neighbouring)
    {
        // courant is defined at each cell center as described by Blazek(2001)
        volScalarField c = sqrt(thermoPtr->Cp()/thermoPtr->Cv()/thermoPtr->psi());
        vectorField magUc = cmptMag(UPtr->internalField()) + c.internalField()*vector::one;
        scalarField lambdaC = magUc & dS();
        rDeltaTDSTPtr->internalField() = lambdaC/(courant * mesh.V());
    }
    else
    {
        // courant is defined as the maximum courant calculated at cell faces (HiSA) 
        volScalarField gamma = thermoPtr->Cp()/thermoPtr->Cv();
        surfaceScalarField amaxSf =
        mag(fvc::interpolate(*UPtr)) +
        fvc::interpolate(sqrt(gamma/thermoPtr->psi()));

        surfaceScalarField dtSf = amaxSf * mesh.deltaCoeffs(); // deltaCoeffs() = 1./max((delta & normal),0.05*mag(delta)) = of6 nonOrthDeltaCoeffs()
        const labelList& owner = mesh.owner();
        const labelList& neighbour=  mesh.neighbour();
        volScalarField& dt = rDeltaTDSTPtr();
        dt *= 0;
        
        forAll(dtSf, faceI)
        {
            dt[owner[faceI]] = max(dt[owner[faceI]], dtSf[faceI]);
            dt[neighbour[faceI]] = max(dt[neighbour[faceI]], dtSf[faceI]);
        }
        
        forAll(dtSf.boundaryField(), patchI)
        {
            const labelList& fc = mesh.boundary()[patchI].faceCells();
            if (mesh.boundary()[patchI].coupled())
            {
                forAll(fc, faceI)
                {
                    dt[fc[faceI]] = max(dt[fc[faceI]], dtSf.boundaryField()[patchI][faceI]);
                }
            }
            else if(mesh.boundary()[patchI].type() == "wall")
            {
                const scalarField deltaP = mesh.deltaCoeffs().boundaryField()[patchI];
                forAll(fc, faceI)
                {
                    scalar pLambda = 0.5*deltaP[faceI]*(sqrt(gamma[fc[faceI]]/thermoPtr->psi()[fc[faceI]]))+mag((*UPtr)[fc[faceI]]);
                    dt[fc[faceI]] = max(dt[fc[faceI]], pLambda);
                }
            }
        }
        dt.internalField() /= courantPtr();
    } 

    Info<<"maxDeltaT "<<(1./gMin(rDeltaTDSTPtr()))<<" minDeltaT "<<(1./gMax(rDeltaTDSTPtr()))<<endl;
}

void timeDST::setTransientDeltaT()
{
    const fvMesh& mesh = UPtr->mesh();
    if(!neighbouring)
    {
        // courant is defined at each cell center as described by Blazek(2001)
        volScalarField c = sqrt(thermoPtr->Cp()/thermoPtr->Cv()/thermoPtr->psi());
        vectorField magUc = cmptMag(UPtr->internalField()) + c.internalField()*vector::one;
        scalarField lambdaC = magUc & dS();
        scalarField deltaTField = courant * mesh.V()/lambdaC;

        setDeltaT
        (
            gMin
            (
                deltaTField
            )
        );
    }
    else
    {
        // courant is defined as the maximum courant calculated at cell faces (HiSA) 
        surfaceScalarField amaxSf =
        mag(fvc::interpolate(*UPtr) & mesh.Sf()) +
        mesh.magSf() * fvc::interpolate(sqrt(thermoPtr->Cp()/thermoPtr->Cv()/thermoPtr->psi()));
        //lambda = fvc::interpolate(sqrt(gamma()/ppsi)) + fvc::interpolate(mag(U_()));    

        surfaceScalarField dtSf = mesh.magSf()/(amaxSf * mesh.deltaCoeffs());        
        
        scalarField deltaTField = dtSf * courant;

        setDeltaT
        (
            gMin
            (
                deltaTField
            )
        );
    }    
}

void timeDST::setCourantDST()
{
    const fvMesh& mesh = UPtr->mesh();
    const volVector5Field& R = mesh.lookupObject<volVector5Field>("R_");
    volScalarField magR = mag(R);
    
    scalar residualR = sqrt(gSumSqr(magR));
    Info<<"residualR = "<<residualR<<endl;

    if(timeIndex_ != 1)
    {
        residualRatio = 
        (
            min
            (
                max
                (
                    residual/residualR,
                    scalar(0.1)
                ),
                scalar(2)
            )
        );
    }
    residual = residualR;    
    
    if (steadyState)
    {
        const volScalarField& rho = mesh.lookupObject<volScalarField>("rho");
        const volVectorField& rhoU = mesh.lookupObject<volVectorField>("rhoU");
        const volScalarField& rhoE = mesh.lookupObject<volScalarField>("rhoE");
        
        // this method is taken from HiSA and offers good performance
        const volScalarField& rhoNew(rho);
        const volScalarField eNew(rhoE/rho-0.5*magSqr(rhoU/rho));
        volScalarField& rhoOld = rhoOldPtr();
        volScalarField& eOld = eOldPtr();
        scalarField factor(mesh.nCells(),scalar(1));
        const labelList& owner = mesh.owner();
        const labelList& neighbour = mesh.neighbour();
        
        forAll(mesh.Sf(), faceI)
        {
            if(rhoNew[owner[faceI]]<0.95*rhoOld[owner[faceI]] || eNew[owner[faceI]]<0.95*eOld[owner[faceI]] || eNew[owner[faceI]]<SMALL)
            {
                factor[owner[faceI]] = min(factor[owner[faceI]],0.5);
                factor[neighbour[faceI]] = min(factor[neighbour[faceI]],0.75);
            }
            
            if(rhoNew[neighbour[faceI]]<0.95*rhoOld[neighbour[faceI]] || eNew[neighbour[faceI]]<0.95*eOld[neighbour[faceI]] || eNew[neighbour[faceI]]<SMALL)
            {
                factor[neighbour[faceI]] = min(factor[neighbour[faceI]],0.5);
                factor[owner[faceI]] = min(factor[owner[faceI]],0.75);
            }
        }

        forAll(mesh.boundary(),patchI)
        {
            if(mesh.boundary()[patchI].coupled())
            {
                scalarField rhoOp = rhoOld.boundaryField()[patchI].patchNeighbourField();
                scalarField rhoNp = rhoNew.boundaryField()[patchI].patchNeighbourField();
                scalarField eOp = eOld.boundaryField()[patchI].patchNeighbourField();
                scalarField eNp = eNew.boundaryField()[patchI].patchNeighbourField();
                const labelUList& fc = mesh.boundary()[patchI].faceCells();
                forAll(fc,faceI)
                {
                    if(rhoNp[faceI]<0.95*rhoOp[faceI] || eNp[faceI]<0.95*eOp[faceI] || eNp[faceI]<SMALL)
                    {
                        factor[fc[faceI]] = min(0.75,factor[fc[faceI]]);
                    }
                }
            }
        }
        rhoOld = rhoNew;
        eOld = eNew;
        
        courantPtr->internalField() *= residualRatio*factor;
        
        courantPtr() = 
        (
            min
            (
                max
                (
                    courantPtr(),
                    minCo
                ),
                maxCo
            )
        );
        
        Info<<"minCourant = "<<min(courantPtr()).value()
            <<" maxCourant = "<<max(courantPtr()).value()
            <<" meanCo = "<<sum(courantPtr()).value()/mesh.globalData().nTotalCells()<<endl;
        setSteadyStateDeltaT();
    }
    else
    {
        courant  = 
        (
            min
            (
                max
                (
                    courant*residualRatio,
                    minCo
                ),
                maxCo
            )
        );
        setTransientDeltaT();
        Info<<"Courant = "<<courant<<", deltaT = "<<deltaTValue()<<endl;
            courantPtr->internalField() = 2.0/3*courant;
    }
    innerIter = 0;
}

bool timeDST::localTimeStepping()
{
    return steadyState;
}
const scalarField& timeDST::rDeltaTDST() const
{
    return rDeltaTDSTPtr->internalField();
}

int timeDST::nIter() const
{
    return innerIter;
}

scalar timeDST::coNum()
{
    return courant;
}

bool timeDST::innerLoop()
{
    innerIter++;
    return 
    (
        innerIter < maxInnerIter
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam
