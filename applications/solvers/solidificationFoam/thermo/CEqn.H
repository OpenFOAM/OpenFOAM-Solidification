{
    // Relative composition for advection/diffusion source terms
    volScalarField CRel("CRel", C2 - C);

    fvScalarMatrix CEqn
    (
        fvm::ddt(rho, C)
      + fvm::div(phi, C)
      - fvm::laplacian(alpha2*rho2*D2, C)
    ==
        fvc::laplacian(alpha2*rho2*D2, CRel)
      - fvc::div(phi, CRel)
    );

    CEqn.relax();
    CEqn.solve();
}
