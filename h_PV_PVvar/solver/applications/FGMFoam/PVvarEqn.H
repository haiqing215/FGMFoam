
if(!laminarBool && kEpsilonBool){
    XiPV = turbulence->epsilon() / (turbulence->k() + smallk) * PVvar; 
}
if(!laminarBool && LESBool){
    scalarField& cellVolumeCells = cellVolumes.primitiveFieldRef();
    cellVolumeCells = mesh.V();
    XiPV = pow( CEpsilon * turbulence->mut() * pow( std::sqrt(2.0) * mag(symm(fvc::grad(U))) , 2) / ( rho * pow(cellVolumes,2.0/3.0) ), 1.0/3.0 ) * PVvar;
}

// Solve the progress variable variance equation
fvScalarMatrix PVvarEqn
(
    fvm::ddt(rho, PVvar)
  + mvConvection->fvmDiv(phi, PVvar)
  - fvc::div((thermo.DPV() + turbulence->mut()/Sct) * fvc::grad(PVvar))
 ==
    2*turbulence->mut()/Sct*magSqr(fvc::grad(PV))
  + rho*XiPV
  + 2*( thermo.PVsourcePV() - PV*thermo.sourcePV())
  + fvOptions(rho, PVvar)
);

PVvarEqn.relax();
fvOptions.constrain(PVvarEqn);
PVvarEqn.solve();

