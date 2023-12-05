// Convection scheme
tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,h)")
    )
);

// Solve the enthalpy equation
fvScalarMatrix hEqn
(
    fvm::ddt(rho, h)
  + mvConvection->fvmDiv(phi, h)
  - fvm::laplacian(thermo.alpha(), h)
 ==
    fvOptions(rho, h)
);

hEqn.relax();
fvOptions.constrain(hEqn);
hEqn.solve();

