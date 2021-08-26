/*-----------------------------------------------------------------------------------------*\
| =========                 |                                                               |
| \\      /  F ield         | rhoDST: Block-coupled solver for compressible flows           |
|  \\    /   O peration     |                                                               |
|   \\  /    A nd           | Copyright (C) 2021 Design and Simulation Tech. Inc. (DSTECH)  |
|    \\/     M anipulation  | Website:  http://www.dstechno.net/                            |
\*-----------------------------------------------------------------------------------------*/
/*
License
    This file is derivative work of rhoDST.

    init-rhoDST is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    init-rhoDST is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFoam.  If not, see <http://www.gnu.org/licenses/>.

Author
-------
    Design and Simulation Technologies Inc. (DSTECH)
    http://dstechno.net/
         
      _____   _____ _______ ______ _____ _    _ 
     |  __ \ / ____|__   __|  ____/ ____| |  | |
     | |  | | (___    | |  | |__ | |    | |__| |
     | |  | |\___ \   | |  |  __|| |    |  __  |
     | |__| |____) |  | |  | |___| |____| |  | |
     |_____/|_____/   |_|  |______\_____|_|  |_|


\*--------------------------------------------------------------------------------------------*/

DST::retrieve(W, rho, 0);
DST::retrieve(W, rhoU, 1);
DST::retrieve(W, rhoE, 4);

dictionary fieldBounds = mesh_.solutionDict().subDict("fieldBounds");

// Pressure bounds
dimensionedScalar pMin("pMin", thermo_.p().dimensions(), SMALL);
dimensionedScalar pMax("pMax", thermo_.p().dimensions(), 1e10);
fieldBounds.lookup(thermo_.p().name()) >> pMin.value() >> pMax.value();

// Temperature bounds
dimensionedScalar TMin("TMin", thermo_.T().dimensions(), SMALL);
dimensionedScalar TMax("TMax", thermo_.T().dimensions(), 1e10);
fieldBounds.lookup(thermo_.T().name()) >> TMin.value() >> TMax.value();

// Density bounds
dimensionedScalar rhoMin("rhoMin", rho.dimensions(), SMALL);
dimensionedScalar rhoMax("rhoMax", rho.dimensions(), 1e10);
fieldBounds.lookup(rho.name()) >> rhoMin.value() >> rhoMax.value();

boundMinMax(rho, rhoMin, rhoMax);

U_ = rhoU/rho;
U_.correctBoundaryConditions();

thermo_.h() = thermo_.Cp()/thermo_.Cv()*(rhoE/rho - 0.5*magSqr(U_)); // h = rhoE/rho - 0.5*magSqr(U) +p/rho; -> unstable 
thermo_.h().correctBoundaryConditions();

dimensionedScalar CpMin = min(thermo_.Cp());
dimensionedScalar CpMax = max(thermo_.Cp());

dimensionedScalar hMin = CpMin*TMin;
dimensionedScalar hMax = CpMax*TMax;

boundMinMax(thermo_.h(), hMin, hMax);

thermo_.correct();
    
//p = rho*(thermo_.Cp() - thermo_.Cv())*T;
thermo_.p().dimensionedInternalField() = rho/thermo_.psi();
thermo_.p().correctBoundaryConditions();

boundMinMax(thermo_.p(), pMin, pMax);

//volScalarField psi2 = (thermo_.Cp() - thermo_.Cv())*T;
rho.boundaryField() = thermo_.p().boundaryField()*thermo_.psi().boundaryField();
rho.correctBoundaryConditions();
rhoU.boundaryField() = rho.boundaryField()*U_.boundaryField();
rhoE.boundaryField() =
    rho.boundaryField()*
    (
        thermo_.h().boundaryField() + 0.5*magSqr(U_.boundaryField())
    ) - thermo_.p().boundaryField();
rhoE.correctBoundaryConditions();
