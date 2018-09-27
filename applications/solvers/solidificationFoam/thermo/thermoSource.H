forAll(mesh.cells(), i)
{
    scalar Thigh = Tsol[i];
    scalar Tlow  = Tsol[i];

    // Aritifical limits for eutectic reaction
    if (C[i] >= Calpha.value())
    {
        Thigh = Teut.value() + Tiso.value();
        Tlow  = Teut.value() - Tiso.value();
    }

    // All liquid
    if (T[i] >= Tliq[i])
    {
        alpha10[i] = 0.0;

        dFdT[i] = 0.0;
        T0[i] = T[i];
    }

    // Primary solidification
    if ((T[i] < Tliq[i]) && (T[i] > Thigh))
    {
        visitedEut[i] = 0;

        dFdT[i] = (Tliq[i] - Tmelt.value()) / (sqr(T[i] - Tmelt.value()))
          / (1.0 - kp.value());


        T0[i] = (Tliq[i] - alpha1[i]*(1.0 - kp.value())*Tmelt.value())
          / (1.0 - alpha1[i]*(1.0 - kp.value()));
    }

    // Possible eutectic solidification (linear regression)
    if ((T[i] <= Thigh) && (T[i] > Tlow))
    {
        if (visitedEut[i] == 0)
        {
            alphaEut[i] = alpha1[i];
            visitedEut[i] = 1;
        }

        dFdT[i] = (alphaEut[i] - 1.0) / (Thigh - Tlow);

        T0[i] = Tlow + (alpha1[i] - 1.0) / dFdT[i];
    }

    // All solid
    if (T[i] <= Tlow)
    {
        alpha10[i] = 1.0;

        dFdT[i] = 0.0;
        T0[i] = T[i];
    }

}

dFdT.correctBoundaryConditions();
T0.boundaryFieldRef() = T.boundaryField();
