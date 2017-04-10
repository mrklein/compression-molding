    Info<< "Reading field T\n" << endl;

    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField tinduction
    (
        IOobject
        (
            "tinduction",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("tinduction", dimless, 0.0)
    );

    volScalarField alpha
    (
        IOobject
        (
            "alpha",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("alpha", dimless, 0.0)
    );

    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );


    Info<< "Reading diffusivity DT\n" << endl;

    dimensionedScalar DT("DT", dimArea/dimTime, transportProperties);
    dimensionedScalar t0("t0", dimTime, transportProperties);
    dimensionedScalar T0("T0", dimTemperature, transportProperties);
    dimensionedScalar k0("k0", dimless/dimTime, transportProperties);
    dimensionedScalar Ea("Ea", dimEnergy/dimMoles, transportProperties);
    dimensionedScalar n("n", dimless, transportProperties);
