{
    // Update mixture properties
    volScalarField Cp(fluid.Cp());
    volScalarField kappa("kappa", fluid.kappa());

    // Semi-implicit handling of eutectic reaction
    volScalarField Sp = rho2*L*dFdT/runTime.deltaT();
    volScalarField Sc = -Sp*T0 - fvc::ddt(rho2*L, alpha2);

    fvScalarMatrix TEqn
    (
        fvm::ddt(rho*Cp, T)
      + fvm::div(phi*fvc::interpolate(Cp), T, "div(phi,T)")
      - fvm::laplacian(kappa, T)
     ==
        fvm::Sp(Sp, T) + Sc
      - L*fvc::div(phi)
    );

    TEqn.relax();
    TEqn.solve();
}
