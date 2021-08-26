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
    const IOobject& iodict
)
:
    IOdictionary(iodict)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

IOdictDST::~IOdictDST()
{}

// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

scalar IOdictDST::beta() const
{
    return timeIntegrationDict().lookupOrDefault<scalar>("beta",1);
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam
