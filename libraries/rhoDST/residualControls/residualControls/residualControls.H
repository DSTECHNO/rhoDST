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
    residualControls

SourceFiles
    residualControls.H

\*---------------------------------------------------------------------------*/
#ifndef residualControls_H
#define residualControls_H

#include "fvCFD.H"
#include "geometricFieldDST.H"


namespace Foam
{

template<class Type>
class geometricFieldDST;

/*---------------------------------------------------------------------------*\
                           Class residualControls Declaration
\*---------------------------------------------------------------------------*/
class residualControls
{
    private:
        
    protected:
        const geometricFieldDST<vector5>& W_;
        const IOdictDST& dict_; 
        DynamicList<scalar> residuals;
        scalar residual;
        int iterationStart;
        int iterationCount;
        
    public:
        //- Constructors
            residualControls(const geometricFieldDST<vector5>& W);

        //- Destructor
            ~residualControls();

        //- Member Functions
            void addResidual(vector5 res);
            int nIter();
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
}

#endif
