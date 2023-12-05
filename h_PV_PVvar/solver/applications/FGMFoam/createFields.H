#include "createRDeltaT.H"

// Thermo model
Info<< "Reading thermophysical properties\n" << endl;
autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
psiReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

// Pressure
volScalarField& p = thermo.p();

// Density
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

// Velocity
Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Access to the control variables
volScalarField& h      = thermo.h();
volScalarField& PV     = thermo.PV();
volScalarField& PVvar  = thermo.PVvar();

// Create interpolation scheme for control variables
multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

fields.add(h);
fields.add(PV);
fields.add(PVvar);

#include "compressibleCreatePhi.H"
pressureControl pressureControl(p, rho, pimple.dict(), false);
mesh.setFluxRequired(p.name());

// Turbulence model
Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

// Set some booleans for the turbulence models 
bool laminarBool  = true;
bool kEpsilonBool = false;
bool LESBool      = false;

// List of all available LES models
wordList LESModels = {"DeardorffDiffStress",
                      "Smagorinsky",
                      "SpalartAllmarasDDES",
                      "SpalartAllmarasDES",
                      "SpalartAllmarasIDDES",
                      "WALE",
                      "dynamicKEqn",
                      "dynamicLagrangian",
                      "kEqn",
                      "kOmegaSSTDDES",
                      "kOmegaSSTDES",
                      "kOmegaSSTIDDES"};

if (turbulence->type() != "Stokes"){ // If not laminar
    laminarBool = false;
    if (turbulence->type() == "kEpsilon"){ // If RANS and kEpsilon
        kEpsilonBool = true;
    }
    else if (std::find(std::begin(LESModels), std::end(LESModels), turbulence->type()) != std::end(LESModels)){ // If LES
        LESBool = true;

        
    }
    else { // Only the RANS kEpsilon and all LES models work with this solver
        FatalErrorInFunction << "\nPlease provide valid turbulence model. Available models are:\n\nRANS:\tkEpsilon\nLES:\tall models" << exit(FatalError);
    }
}

#include "createDpdt.H"

#include "createK.H"

#include "createMRF.H"

const dimensionedScalar rhoMax("rhoMax", dimDensity, GREAT, pimple.dict());
const dimensionedScalar rhoMin("rhoMin", dimDensity, Zero, pimple.dict());

#include "createFvOptions.H"

// Combustion model
Info<< "Creating combustion model\n" << endl;
autoPtr<CombustionModel<psiReactionThermo>> combustion
(
    CombustionModel<psiReactionThermo>::New(thermo, turbulence())
);

const word combModelName(combustion->type());
Info<< "Combustion Model Used: " << combModelName << endl;

// Density
rho = thermo.rho();

// Turbulent scalar dissipation of the progress variable
volScalarField XiPV
(
    IOobject
    (
        "XiPV",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimensionSet(0,0,-1,0,0,0,0), Zero)
);

// Turbulent Schmidt number
volScalarField Sct
(
    IOobject
    (
        "Sct",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, 0.7)
);

// LES turbulent scalar dissipation model parameter
volScalarField CEpsilon
(
    IOobject
    (
        "CEpsilon",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, 1.8)
);

// Thermal diffusivity
volScalarField alpha
(
    IOobject
    (
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.alpha()
);

// Dynamic viscosity
volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.mu()
);

// Arbritrary small k value to prevent division by zero
volScalarField smallk
(
    IOobject
    (
        "smallk",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimensionSet(0,2,-2,0,0,0,0), 1E-10)
);

// Field with cell volumes. This is needed to ensure type coherence in the Xi computation
volScalarField cellVolumes
(
    IOobject
    (
        "cellVolumes",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimensionSet(0,3,0,0,0,0,0), 1E-10)
);

// Field with change in PV per cell. This is needed for the adaptive mesh refinement.
volScalarField deltaPV
(
    IOobject
    (
        "deltaPV",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimensionSet(0,0,0,0,0,0,0), 1E-10)
);
