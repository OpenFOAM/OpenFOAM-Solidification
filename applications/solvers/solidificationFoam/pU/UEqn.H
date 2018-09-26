MRF.correctBoundaryVelocity(U);

// Calculate the drag force
volScalarField KInv
(
    (180.0/sqr(DAS))*pow(min(alpha1, 0.95), 2.0)/pow(max(alpha2, 0.05), 3.0)
);

volScalarField Kd("Kd", rho*mu2/rho2*KInv);

// Update the kinematic density
rhok = rho2*(1.0 - betaT*(T - TRef) - betaC*(C2 - CRef));


fvVectorMatrix UEqn
(
    fvm::ddt(rho, U) + fvm::div(phi, U)
  + MRF.DDt(rho, U)
  - fvm::laplacian(rho*mu2/rho2, U)
 ==
    fvOptions(rho, U)
  - fvm::Sp(Kd, U)
);
UEqn.relax();
fvOptions.constrain(UEqn);


if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
      ==
        fvc::reconstruct
        (
            (
              - ghf*fvc::snGrad(rhok)
              - fvc::snGrad(p_rgh)
            )*mesh.magSf()
        )
    );

    fvOptions.correct(U);
}
