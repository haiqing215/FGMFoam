/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "fluidThermo.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(fluidThermo, 0);
    defineRunTimeSelectionTable(fluidThermo, fvMesh);
    defineRunTimeSelectionTable(fluidThermo, fvMeshDictPhase);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidThermo::fluidThermo(const fvMesh& mesh, const word& phaseName)
:
    basicThermo(mesh, phaseName),

    h_
    (
        IOobject
        (
            phasePropertyName("h"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
    mesh
    ),

    PV_
    (
        IOobject
        (
            phasePropertyName("PV"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
    mesh
    ),

    PVvar_
    (
        IOobject
        (
            phasePropertyName("PVvar"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
    mesh
    ),

    sourcePV_
    (
        IOobject
        (
            phasePropertyName("sourcePV"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimDensity/dimTime, 0)
    ),

    DPV_
    (
        IOobject
        (
            phasePropertyName("DPV"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimMass/dimLength/dimTime, 0)
    ),

    PVsourcePV_
    (
        IOobject
        (
            phasePropertyName("PVsourcePV"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimDensity/dimTime, 0)
    )

{}



Foam::fluidThermo::fluidThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    basicThermo(mesh, dict, phaseName),

    h_
    (
        IOobject
        (
            phasePropertyName("h"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
    mesh
    ),

    PV_
    (
        IOobject
        (
            phasePropertyName("PV"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
    mesh
    ),

    PVvar_
    (
        IOobject
        (
            phasePropertyName("PVvar"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
    mesh
    ),

    sourcePV_
    (
        IOobject
        (
            phasePropertyName("sourcePV"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimDensity/dimTime, 0)
    ),

    DPV_
    (
        IOobject
        (
            phasePropertyName("DPV"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimMass/dimLength/dimTime, 0)
    ),

    PVsourcePV_
    (
        IOobject
        (
            phasePropertyName("PVsourcePV"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimDensity/dimTime, 0)
    )

{}


Foam::fluidThermo::fluidThermo
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictionaryName
)
:
    basicThermo(mesh, phaseName, dictionaryName),

    h_
    (
        IOobject
        (
            phasePropertyName("h"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
    mesh
    ),

    PV_
    (
        IOobject
        (
            phasePropertyName("PV"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
    mesh
    ),

    PVvar_
    (
        IOobject
        (
            phasePropertyName("PVvar"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
    mesh
    ),

    sourcePV_
    (
        IOobject
        (
            phasePropertyName("sourcePV"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimDensity/dimTime, 0)
    ),

    DPV_
    (
        IOobject
        (
            phasePropertyName("DPV"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimMass/dimLength/dimTime, 0)
    ),

    PVsourcePV_
    (
        IOobject
        (
            phasePropertyName("PVsourcePV"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
    dimensionedScalar(dimDensity/dimTime, 0)
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fluidThermo> Foam::fluidThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<fluidThermo>(mesh, phaseName);
}


Foam::autoPtr<Foam::fluidThermo> Foam::fluidThermo::New
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictName
)
{
    return basicThermo::New<fluidThermo>(mesh, phaseName, dictName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidThermo::~fluidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::fluidThermo::nu() const
{
    return mu()/rho();
}


Foam::tmp<Foam::scalarField> Foam::fluidThermo::nu(const label patchi) const
{
    return mu(patchi)/rho(patchi);
}

Foam::volScalarField& Foam::fluidThermo::h()
{
    return h_;
}

const Foam::volScalarField& Foam::fluidThermo::h() const
{
    return h_;
}

Foam::volScalarField& Foam::fluidThermo::PV()
{
    return PV_;
}

const Foam::volScalarField& Foam::fluidThermo::PV() const
{
    return PV_;
}

Foam::volScalarField& Foam::fluidThermo::PVvar()
{
    return PVvar_;
}

const Foam::volScalarField& Foam::fluidThermo::PVvar() const
{
    return PVvar_;
}

Foam::volScalarField& Foam::fluidThermo::sourcePV()
{
    return sourcePV_;
}

const Foam::volScalarField& Foam::fluidThermo::sourcePV() const
{
    return sourcePV_;
}

Foam::volScalarField& Foam::fluidThermo::DPV()
{
    return DPV_;
}

const Foam::volScalarField& Foam::fluidThermo::DPV() const
{
    return DPV_;
}

Foam::volScalarField& Foam::fluidThermo::PVsourcePV()
{
    return PVsourcePV_;
}

const Foam::volScalarField& Foam::fluidThermo::PVsourcePV() const
{
    return PVsourcePV_;
}

// ************************************************************************* //
