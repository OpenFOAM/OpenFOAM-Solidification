for (int tCorr=0; tCorr<nThermoCorr; tCorr++)
{
    // Store previous values for linearized update
    volScalarField alpha10(alpha1);
    volScalarField T0(T);

    // Update the freezing range
    volScalarField Tliq(min(max(Tmelt + mliq*C, Teut), Tmelt));
    volScalarField Tsol(min(max(Tmelt + msol*C, Teut), Tmelt));

    #include "thermo/thermoSource.H"
    #include "thermo/TEqn.H"
    #include "thermo/thermoUpdate.H"
    #include "thermo/CEqn.H"

    // Check for phase fraction convergence
    scalarField residualField
    (
        mag(alpha1.internalField() - alpha10.internalField())
    );

    scalar residual(gMax(residualField));

    Info<< "THERMO: iteration " << (tCorr + 1)
        << " residual: " << residual  << endl;

    if (residual < 1e-4)
    {
        break;
    }
}
