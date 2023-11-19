/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "baseFGM.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::baseFGM<ReactionThermo>::baseFGM
(
    const word& modelType,
    ReactionThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    ThermoCombustion<ReactionThermo>(modelType, thermo, turb),
    PV_         (this->thermo().PV()),
    PVvar_      (this->thermo().PVvar()),
    sourcePV_   (this->thermo().sourcePV()),
    DPV_        (this->thermo().DPV()),
    PVsourcePV_ (this->thermo().PVsourcePV()),
    T_          (this->thermo().T()),
    mu_         (const_cast<volScalarField&>(this->thermo().mu()())),
    alpha_      (const_cast<volScalarField&>(this->thermo().alpha())),
    psi_        (const_cast<volScalarField&>(this->thermo().psi())),
      
    fgm_(readFGM("constant/lookUp/database.fgm"))
{
    // Initialise control and total variables as null pointer (empty)
    controlVariables_ = nullptr;
    variables_        = nullptr;

    // Release memory if arrays already exist
    if (controlVariables_) {
        delete[] controlVariables_;
    }
    if (controlVariables_) {
        delete[] controlVariables_;
    }

    // Initialize arrays with sizes obtained from the FGM table
    controlVariables_ = new double[fgm_->Ncv];
    variables_        = new double[fgm_->Nvar];

    // Update the field variables with the table variables
    update();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::baseFGM<ReactionThermo>::~baseFGM()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ReactionThermo>
void Foam::combustionModels::baseFGM<ReactionThermo>::update()
{ 
    const clockTime clockTime_= clockTime();
    clockTime_.timeIncrement();
    
    scalar PVmin = 0.0;
    scalar PVmax = 7.801199268062700e-02;  // HARDCODED
   
    // Initialise cp and lambda at arbitrary small values. These will be overwritten by the table data 
    scalar cp     = 1.0e-10;
    scalar lambda = 1.0e-10;

    // Create reference to the field cell center values
    scalarField& PVCells         = this->PV_.primitiveFieldRef();
    scalarField& PVvarCells      = this->PVvar_.primitiveFieldRef();

    scalarField& sourcePVCells   = this->sourcePV_.primitiveFieldRef();
    scalarField& DPVCells        = this->DPV_.primitiveFieldRef();
    scalarField& PVsourcePVCells = this->PVsourcePV_.primitiveFieldRef();
    scalarField& TCells          = this->T_.primitiveFieldRef();
    scalarField& muCells         = this->mu_.primitiveFieldRef();
    scalarField& alphaCells      = this->alpha_.primitiveFieldRef();
    scalarField& psiCells        = this->psi_.primitiveFieldRef();

    // Loop over all cells
    forAll(PVCells, celli)
    {
        // Scale the control variables
        controlVariables_[0] = (PVCells[celli] - PVmin)/(PVmax - PVmin);
        controlVariables_[1] = ((PVvarCells[celli] + pow(PVCells[celli],2) - pow(PVmin,2) - 2*controlVariables_[0]*(PVmin*PVmax - pow(PVmin,2)))/(pow((PVmax-PVmin),2))) - pow(controlVariables_[0],2);
        
        // Edge treatment
        if (controlVariables_[0] > 1.0) {
            controlVariables_[0] = 1.0;
        }
        if (controlVariables_[0] < 0.0) {
            controlVariables_[0] = 0.0;
        }
        if (controlVariables_[1] > 0.25) {
            controlVariables_[1] = 0.25;
        }
        if (controlVariables_[1] < 0.0) {
            controlVariables_[1] = 0.0;
        }

        lookupFGM_2D(fgm_,controlVariables_,variables_);

        for (int i = fgm_->Ncv; i < fgm_->Nvar; i++)
        {
            if (std::strcmp(fgm_->varname[i], "SOURCE_CV1") == 0){
                sourcePVCells[celli] = variables_[i];
            }
            if (std::strcmp(fgm_->varname[i], "YOMEGAY") == 0){
                PVsourcePVCells[celli] = variables_[i];
            }
            if (std::strcmp(fgm_->varname[i], "DIFF_CV1") == 0){
                DPVCells[celli] = variables_[i];
            }
            if (std::strcmp(fgm_->varname[i], "TEMPERATURE") == 0){
                TCells[celli] = variables_[i];
            }
            if (std::strcmp(fgm_->varname[i], "VISCOSITY") == 0){
                muCells[celli] = variables_[i];
            }
            if (std::strcmp(fgm_->varname[i], "CP") == 0){
                cp = variables_[i];
            }
            if (std::strcmp(fgm_->varname[i], "CONDUCTIVITY") == 0){
                lambda = variables_[i];
            }
            if (std::strcmp(fgm_->varname[i], "DENSITY") == 0){
                psiCells[celli] = variables_[i]/101325.0;
            }
        }

        // If cp and lambda are in the variable list, update alpha
        if (cp != 0 && lambda != 0){
            alphaCells[celli] = lambda/cp;
        }
        else {
            Info << "WARNING: Cp and/or lambda do not have a value, alpha does not have a value" << endl;
        }

        // Set sourceterm at zero when network limits are exceeded
        if (PVCells[celli] <= PVmin || PVCells[celli] >= PVmax)
        {
            sourcePVCells[celli]    = 0.0;
            PVsourcePVCells[celli]  = 0.0;
        }

        // Source term is strictly positive
        if (sourcePVCells[celli] < 0.0)
        {
            sourcePVCells[celli]   = 0.0;
        }
        if (PVsourcePVCells[celli] < 0.0)
        {
            PVsourcePVCells[celli]   = 0.0;
        }
    }

    forAll(T_.boundaryFieldRef(), patchi)
    {
        fvPatchScalarField& pPV         = this->PV_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pPVvar      = this->PVvar_.boundaryFieldRef()[patchi];
        
        fvPatchScalarField& psourcePV   = this->sourcePV_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pDPV        = this->DPV_.boundaryFieldRef()[patchi]; 
        fvPatchScalarField& pPVsourcePV = this->PVsourcePV_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pT          = this->T_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pmu         = this->mu_.boundaryFieldRef()[patchi];
        fvPatchScalarField& palpha      = this->alpha_.boundaryFieldRef()[patchi];
        fvPatchScalarField& ppsi        = this->psi_.boundaryFieldRef()[patchi];

        forAll(pT, facei)
        {
            controlVariables_[0] = (pPV[facei] - PVmin)/(PVmax - PVmin);
            controlVariables_[1] = ((pPVvar[facei] + pow(pPV[facei],2) - pow(PVmin,2) - 2*controlVariables_[0]*(PVmin*PVmax - pow(PVmin,2)))/(pow((PVmax-PVmin),2))) - pow(controlVariables_[0],2);
 
            // Edge treatment
            if (controlVariables_[0] > 1.0) {
                controlVariables_[0] = 1.0;
            }
            if (controlVariables_[0] < 0.0) {
                controlVariables_[0] = 0.0;
            }
            if (controlVariables_[1] > 0.25) {
                controlVariables_[1] = 0.25;
            }
            if (controlVariables_[1] < 0.0) {
                controlVariables_[1] = 0.0;
            }

            lookupFGM_2D(fgm_,controlVariables_,variables_);
            
            for (int i = fgm_->Ncv; i < fgm_->Nvar; i++)
            {
                if (std::strcmp(fgm_->varname[i], "DIFF_CV1") == 0){
                    pDPV[facei] = variables_[i];
                }
                if (std::strcmp(fgm_->varname[i], "TEMPERATURE") == 0){
                    if (!pT.fixesValue()){
                        pT[facei] = variables_[i];
                    }
                }
                if (std::strcmp(fgm_->varname[i], "VISCOSITY") == 0){
                    pmu[facei] = variables_[i];
                }
                if (std::strcmp(fgm_->varname[i], "CP") == 0){
                    cp = variables_[i];
                }
                if (std::strcmp(fgm_->varname[i], "CONDUCTIVITY") == 0){
                    lambda = variables_[i];
                }
                if (std::strcmp(fgm_->varname[i], "DENSITY") == 0){
                    ppsi[facei] = variables_[i]/101325.0;
                }
            }
            
            // No reaction at the boundary faces
            psourcePV[facei]   = 0.0;
            pPVsourcePV[facei] = 0.0;
            
            // If cp and lambda are in the variable list, update alpha
            if (cp != 0 && lambda != 0){
                palpha[facei] = lambda/cp;
            }
            else {
                Info << "WARNING: Cp and/or lambda do not have a value, alpha does not have a value" << endl;
            }
        }
    }
    
    Info << "Parameter update time: (" << clockTime_.timeIncrement() << " s)" << endl;
}


template<class ReactionThermo>
void Foam::combustionModels::baseFGM<ReactionThermo>::correct()
{
    update();
}


template<class ReactionThermo>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::baseFGM<ReactionThermo>::R
(
    volScalarField& Y
) const
{
    tmp<fvScalarMatrix> tSu
    (
        new fvScalarMatrix(Y, dimMass/dimTime)
    );

    return tSu;
}


template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::baseFGM<ReactionThermo>::Qdot() const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            this->thermo().phasePropertyName(typeName + ":Qdot"),
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh(),
        dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
    );
}

// ************************************************************************* //
