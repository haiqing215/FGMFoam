/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::combustionModels::baseFGM

Group
    grpCombustionModels

Description
    Dummy combustion model for 'no combustion'

SourceFiles
    baseFGM.C

\*---------------------------------------------------------------------------*/

#ifndef baseFGM_H
#define baseFGM_H

#include "ThermoCombustion.H"
#include "fvOptions.H"
#include "hashedWordList.H"
#include "clockTime.H"

extern "C"
{
    #include "FGMlib.h"
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                        Class baseFGM Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo>
class baseFGM
:
    public ThermoCombustion<ReactionThermo>
{

protected:

    // Protected data

        // Control variable related variables

            //- Enthalpy
            volScalarField& h_;
            
            //- Progress variable
            volScalarField& PV_;
            
            //- Progress variable
            volScalarField& PVvar_;
            
            //- Progress variable source term
            volScalarField& sourcePV_;

            //- Progress variable preferential diffusion molecular diffusion coefficient
            volScalarField& DPV_;

            //- Progress variable multiplied by its source term
            volScalarField& PVsourcePV_;

        // Thermo variables

            //- Temperature
            volScalarField& T_;

            //- Dynamic viscosity
            volScalarField& mu_;

            //- Thermal diffusivity
            volScalarField& alpha_;

            //- Compressibility
            volScalarField& psi_;

        // Scaled control variables for look up

            // - Scaled progress variable
            //volScalarField CV1_;

            // - Scaled progress variable variance
            //volScalarField CV1var_;
 
        // Look-up

            //- Create pointer to the FGM look-up routine
            FGM* fgm_;

            //- Control variables
            double* controlVariables_;

            //- Total variables
            double* variables_;

    // Protected Member Functions
        
        //- Update the thermophysical parameters
        virtual void update();

public:

    //- Runtime type information
    TypeName("baseFGM");


    // Constructors

        //- Construct from components
        baseFGM
        (
            const word& modelType,
            ReactionThermo& thermo,
            const compressibleTurbulenceModel& turb,
            const word& combustionProperties
        );


    //- Destructor
    virtual ~baseFGM();


    // Member Functions

        //- Correct the variables by solving equations and updating thermophysical parameters
        virtual void correct(); 
        
        //- Fuel consumption rate matrix
        virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

        //- Heat release rate [kg/m/s3]
        virtual tmp<volScalarField> Qdot() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "baseFGM.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
