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

    rhoDST is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    init-rhoDST is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFoam.  If not, see <http://www.gnu.org/licenses/>.

Application
-----------
    rhoDST

Description
-----------
    A fully-coupled solver for the solution of high speed compressible flows over 
    aeronautical flows. 

Author:
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

#include "fvCFD.H"
#include "Vector5typedefs.H"
#include "fvBlockMatrix.H"
#include "basicPsiThermo.H"
#include "turbulenceModel.H"
#include "includeDST.H"
#include "IOReferencer.H"
#include "geometricFieldDST.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * *//

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTimeDST.H"
#   include "createMesh.H"
#   include "createDSTDict.H"
#   include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * *//

    Info<< "\nStarting time loop for rhoDST\n" << endl;   

    while (runTime.run())
    {
        runTime++;
        
        rhoDST.iterate
        (
            rho,
            rhoU,
            rhoE
        );
        
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }    

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
