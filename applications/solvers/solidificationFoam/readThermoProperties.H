IOdictionary thermoProperties
(
    IOobject
    (
        "thermoProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar Tmelt
(
    "Tmelt",
    dimTemperature,
    thermoProperties.lookup("Tmelt")
);

dimensionedScalar Teut
(
    "Teut",
    dimTemperature,
    thermoProperties.lookup("Teut")
);

dimensionedScalar Tiso
(
    "Tiso",
    dimTemperature,
    thermoProperties.lookup("Tiso")
);

dimensionedScalar Calpha
(
    "Calpha",
    dimless,
    thermoProperties.lookup("Calpha")
);

dimensionedScalar Ceut
(
    "Ceut",
    dimless,
    thermoProperties.lookup("Ceut")
);

dimensionedScalar L
(
    "L",
    dimEnergy/dimMass,
    thermoProperties.lookup("L")
);

// solidus and liquidus slopes
dimensionedScalar msol = (Teut - Tmelt)/Calpha;
dimensionedScalar mliq = (Teut - Tmelt)/Ceut;

dimensionedScalar kp = Calpha/Ceut;

const scalar nThermoCorr
(
    pimple.dict().lookupOrDefault<scalar>("nThermoCorr", 100)
);
