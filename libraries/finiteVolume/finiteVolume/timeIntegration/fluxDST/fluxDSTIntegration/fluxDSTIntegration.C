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

namespace Foam
{

namespace fvDST
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class Type>
fluxDSTIntegration<Type>::fluxDSTIntegration()
{} // End of constructor

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
template<class Type>
fluxDSTIntegration<Type>::~fluxDSTIntegration()
{} // End of destructor

// * * * * * * * * * * * * * Member Function  * * * * * * * * * * * * * * * //

template<class Type>
tmp<fluxDSTIntegration<Type>> fluxDSTIntegration<Type>::New
(
    const IOdictDST& dict
)
{
    word F = "fluxDST";

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_ -> find(F);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "timeIntegration::New()"
        )   << "timeIntegration: " << F
            << " is not a valid choice. " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return cstrIter()();
}

}

} // End namespace Foam

