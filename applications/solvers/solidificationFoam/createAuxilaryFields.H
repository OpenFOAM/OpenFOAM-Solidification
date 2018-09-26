// Liquidus temperature
volScalarField Tliq
(
    IOobject
    (
        "Tliq",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Tliq", dimTemperature, 0)
);

// Solidus temperature
volScalarField Tsol
(
    IOobject
    (
        "Tsol",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Tsol", dimTemperature, 0)
);

// Phase compositions
volScalarField C1
(
    IOobject
    (
        IOobject::groupName("C", "solid"),
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    C
);

volScalarField C2
(
    IOobject
    (
        IOobject::groupName("C", "liquid"),
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    C
);


// Special handling of eutectic reaction
volScalarField::Internal visitedEut
(
    IOobject
    (
        "visitedEut",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("visitedEut", dimless, 0)
);

volScalarField::Internal alphaEut
(
    IOobject
    (
        "alphaEut",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("alphaEut", dimless, 0.0)
);

volScalarField dFdT
(
    IOobject
    (
        "dFdT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("dFdT", dimless/dimTemperature, 0.0)
);
