// Solve the progress variable equation
fvScalarMatrix PVEqn
(
    fvm::ddt(rho, PV)
  + mvConvection->fvmDiv(phi, PV)
  - fvc::div((thermo.DPV() + turbulence->mut()/Sct) * fvc::grad(PV))
 ==
    thermo.sourcePV()
  + fvOptions(rho, PV)
);

PVEqn.relax();
fvOptions.constrain(PVEqn);
PVEqn.solve();

