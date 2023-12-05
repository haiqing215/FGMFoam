/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::lookUpEnthalpyFvPatchScalarField

Group
    grpInletBoundaryConditions grpOutletBoundaryConditions

Description
    This boundary condition sets the enthalpy boundary values to the values 
    obtained from a neural network.

Usage

SourceFiles
    lookUpEnthalpyFvPatchScalarField.C

See also
    Foam::fixedValueFvPatchField

\*---------------------------------------------------------------------------*/

#ifndef lookUpEnthalpyFvPatchScalarField_H
#define lookUpEnthalpyFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

extern "C"
{
    #include "FGMlib.h"
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class lookUpEnthalpyFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class lookUpEnthalpyFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Create pointer to the FGM look-up routine
        FGM* fgm_;

        //- Control variables
        double* controlVariables_;

        //- Total variables
        double* variables_;

        //- List of control variable values
        //List<fvPatchScalarField> controlVariables_;
        //const List<volScalarField*> controlVariables;
        
        //- Total temperature
        scalarField h0_;
        
        //- Total temperature
        //scalarField T0_;

public:

    //- Runtime type information
    TypeName("lookUpEnthalpy");


    // Constructors

        //- Construct from patch and internal field
        lookUpEnthalpyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        lookUpEnthalpyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given lookUpEnthalpyFvPatchScalarField
        //  onto a new patch
        lookUpEnthalpyFvPatchScalarField
        (
            const lookUpEnthalpyFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        lookUpEnthalpyFvPatchScalarField
        (
            const lookUpEnthalpyFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new lookUpEnthalpyFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        lookUpEnthalpyFvPatchScalarField
        (
            const lookUpEnthalpyFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new lookUpEnthalpyFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the enthalpy
            const scalarField& h0() const
            {
                return h0_;
            }

            //- Return reference to the enthalpy to allow adjustment
            scalarField& h0()
            {
                return h0_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
