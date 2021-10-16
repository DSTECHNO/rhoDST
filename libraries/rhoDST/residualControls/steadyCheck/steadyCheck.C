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

#include "steadyCheck.H"
namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

steadyCheck::steadyCheck
(
    const geometricFieldDST<vector5>& W
)
:
    residualControls(W)
{
    magConvergedSlope = dict_.magConvergedSlope();
    minResudual = dict_.minResudual();
    magScaledSlope = 1;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

steadyCheck::~steadyCheck()
{}

// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //


bool steadyCheck::converged()
{
    if(iterationCount>=iterationStart)
    {
        scalar sumXi = 0;
        scalar sumYi = 0;
        scalar sumXiSqr = 0;
        scalar sumYiXi = 0;
        int N = 0.2*residuals.size();
        scalar maxY = log(max(residuals));
        scalar minY = log(min(residuals));
        scalar deltaY = maxY - minY;
        for(int i = 0; i<N; i++)
        {
            int xi = residuals.size()-1-i;
            scalar yi = log(residuals[xi]);
            sumXi += xi;// TODO: use algebraic expressions if possible
            sumXiSqr += pow(xi,2);
            sumYi += yi;
            sumYiXi += yi*xi;
        }
        magScaledSlope = (residuals.size()-1)/deltaY*mag((sumYiXi - sumXi*sumYi/N)/(sumXiSqr-pow(sumXi,2)/N));
    }
    return (iterationCount > dict_.maxInnerIter() || (residual<minResudual || magScaledSlope < magConvergedSlope));
}

bool steadyCheck::converged() const // called in TimeDST class
{
    return (iterationCount > dict_.maxInnerIter() || (residual<minResudual || magScaledSlope < magConvergedSlope));
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam
