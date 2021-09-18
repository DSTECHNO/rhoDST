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
    rhoDST::IOdictDST

SourceFiles
    IOdictDST.C

\*---------------------------------------------------------------------------*/
#ifndef IOdictDST_H
#define IOdictDST_H

#include "IOdictionary.H"
#include "fvCFD.H"


namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class IOdictDST Declaration
\*---------------------------------------------------------------------------*/
class IOdictDST
:
    public IOdictionary
{
    private:
        word forcedSchemeName;
        ITstream* ddtSchemeRho_;
        ITstream* ddtSchemeRhoU_;
        ITstream* ddtSchemeRhoE_;
        const fvMesh& mesh_;
        word forcesScheme;
    protected:

    public:
        //- Constructors
            IOdictDST
            (
                const IOobject& iodict,
                const fvMesh& mesh
            );

        //- Destructor
            ~IOdictDST();

        //- Member Functions
            void setDdtScheme(word ddtScheme);
            bool forcedScheme();
            bool steadyState() const;
            ITstream& ddtSchemeRho();
            ITstream& ddtSchemeRhoU();
            ITstream& ddtSchemeRhoE();
            scalar beta() const;
            bool innerCorrectTurbulence() const;
            bool dualTime() const;
            bool adjustTimeStep() const;
            word dependentVariable() const;
            word fluxScheme() const;
            word Jacobian() const;
            const dictionary& timeIntegrationDict() const;
            const dictionary& timeControlsDict() const;
            const dictionary& residualControlsDict() const;
            word timeIntegration() const;
            bool neighbouring() const;
            scalar initialCo() const;
            scalar minCo() const;
            scalar maxCo() const;
            int iterationStart() const;
            scalar magConvergedSlope() const;
            scalar minResudual() const;
            int maxInnerIter() const;
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
}

#endif
