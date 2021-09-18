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

#include "IOdictDST.H"
namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

IOdictDST::IOdictDST
(
    const IOobject& iodict,
    const fvMesh& mesh
)
:
    IOdictionary(iodict),
    mesh_(mesh)
{
    forcedSchemeName = "default";
    ddtSchemeRho_ = new ITstream("ddtRho",List<token>(1,token(word(forcedSchemeName),1)));
    ddtSchemeRhoU_ = new ITstream("ddtRhoU",List<token>(1,token(word(forcedSchemeName),1)));
    ddtSchemeRhoE_ = new ITstream("ddtRhoE",List<token>(1,token(word(forcedSchemeName),1)));
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

IOdictDST::~IOdictDST()
{
    if(forcedScheme())
    {
        delete ddtSchemeRho_;
        delete ddtSchemeRhoU_;
        delete ddtSchemeRhoE_;
    }
}

// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //


void IOdictDST::setDdtScheme(word ddtScheme)
{
    forcedSchemeName = ddtScheme;
}

bool IOdictDST::forcedScheme()
{
    return forcedSchemeName != "default";
}

bool IOdictDST::steadyState() const
{
    return forcedSchemeName == "default" ? word(mesh_.schemesDict().ddtScheme("rho")) == "steadyStateDST" : forcedSchemeName == "steadyStateDST";
}

ITstream& IOdictDST::ddtSchemeRho()
{
    if(forcedScheme())
    {
        *ddtSchemeRho_ = ITstream("ddtRho",List<token>(1,token(word(forcedSchemeName),1)));
        return *ddtSchemeRho_;
    }

    return mesh_.schemesDict().ddtScheme("rho");
}

ITstream& IOdictDST::ddtSchemeRhoU()
{
    if(forcedScheme())
    {
        *ddtSchemeRhoU_ = ITstream("ddtRhoU",List<token>(1,token(word(forcedSchemeName),1)));
        return *ddtSchemeRhoU_;
    }

    return mesh_.schemesDict().ddtScheme("rhoU");
}

ITstream& IOdictDST::ddtSchemeRhoE()
{
    if(forcedScheme())
    {
        *ddtSchemeRhoE_ = ITstream("ddtRhoE",List<token>(1,token(word(forcedSchemeName),1)));
        return *ddtSchemeRhoE_;
    }

    return mesh_.schemesDict().ddtScheme("rhoE");
}

scalar IOdictDST::beta() const
{
    return timeIntegrationDict().lookupOrDefault<scalar>("beta",1);
}

bool IOdictDST::innerCorrectTurbulence() const
{
    return timeIntegrationDict().lookupOrDefault<bool>("innerCorrectTurbulence",false);
}

bool IOdictDST::adjustTimeStep() const
{
    return timeControlsDict().lookupOrDefault<bool>("adjustTimeStep",true);
}

bool IOdictDST::dualTime() const
{
    return timeIntegration() == "Dual";
}

word IOdictDST::dependentVariable() const
{
    return this->lookupOrDefault<word>("dependentVariable","W");
}

word IOdictDST::fluxScheme() const
{
    return this->lookup("fluxScheme");
}

word IOdictDST::Jacobian() const
{
    return this->lookupOrDefault<word>("Jacobian","Roe");
}

const dictionary& IOdictDST::timeIntegrationDict() const
{
    return this->subDict("timeIntegration");
}

const dictionary& IOdictDST::timeControlsDict() const
{
    return this->subDict("timeControls");
}

const dictionary& IOdictDST::residualControlsDict() const
{
    return this->subDict("residualControls");
}

word IOdictDST::timeIntegration() const
{
    return timeIntegrationDict().lookupOrDefault<word>("timeIntegration","Single");
}

bool IOdictDST::neighbouring() const
{
    return timeControlsDict().lookupOrDefault<bool>("neighbouring",false);
}

scalar IOdictDST::initialCo() const
{
    return timeControlsDict().lookupOrDefault<scalar>("initialCo",1);
}

scalar IOdictDST::minCo() const
{
    return timeControlsDict().lookupOrDefault<scalar>("minCo",.1);
}

scalar IOdictDST::maxCo() const
{
    return timeControlsDict().lookupOrDefault<scalar>("maxCo",50);
}

int IOdictDST::iterationStart() const
{
    int val = 0;
    if(steadyState())
        val = 500;
    else
        val = 10;
    return residualControlsDict().lookupOrDefault<int>("iterationStart",val);
}

scalar IOdictDST::magConvergedSlope() const
{
    return residualControlsDict().lookupOrDefault<scalar>("magConvergedSlope",.1);
}

scalar IOdictDST::minResudual() const
{
    return residualControlsDict().lookupOrDefault<scalar>("minResudual",1e-5);
}

int IOdictDST::maxInnerIter() const
{
    return timeIntegrationDict().lookupOrDefault<int>("maxInnerIter",20);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam
