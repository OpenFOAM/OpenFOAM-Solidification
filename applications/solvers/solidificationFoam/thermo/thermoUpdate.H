forAll(mesh.cells(), i)
{
    // Update the phase fraction due to heat/solute transport
    alpha1[i] = alpha10[i] + dFdT[i]*(T[i] - T0[i]);

    alpha1[i] = min(max(alpha1[i], 0.0), 1.0);

    // Update the phase compositions
    C2[i] = C[i] / (1.0 + alpha1[i]*(kp.value() - 1.0));

    C2[i] = min(C2[i], Ceut.value());
    C1[i] = kp.value()*C2[i];
}

alpha1.correctBoundaryConditions();
alpha2 = scalar(1) - alpha1;

C1.correctBoundaryConditions();
C2.correctBoundaryConditions();
