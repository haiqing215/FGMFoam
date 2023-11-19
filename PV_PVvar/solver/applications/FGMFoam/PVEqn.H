// Convection scheme
tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,PV)")
    )
);

// Solve the progress variable equation
fvScalarMatrix PVEqn
(
    fvm::ddt(rho, PV)
  + mvConvection->fvmDiv(phi, PV)
  - fvc::div((DPV + turbulence->mut()/Sct) * fvc::grad(PV))
 ==
    sourcePV
  + fvOptions(rho, PV)
);

PVEqn.relax();
fvOptions.constrain(PVEqn);
PVEqn.solve();

