// This file is base of OpenFOAM. Simple modification of laplacianFoam.

#include "fvCFD.H"
#include "simpleControl.H"
#include "constants.H"


int main(int argc, char *argv[])
{
    using constant::physicoChemical::R;

    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        Info<< "\nCalculating temperature distribution\n" << endl;

        while (simple.correctNonOrthogonal())
        {
            solve
            (
                fvm::ddt(T) - fvm::laplacian(DT, T)
            );
        }

        volScalarField ti(t0*Foam::exp(T0/T));
        tinduction += runTime.deltaT()/ti;

        volScalarField xi(pos(1 - tinduction));
        volScalarField kprime(Foam::pow(k0*Foam::exp(-Ea/R/T), 1.0/n));

        Info<< "\nCalculating cure rate distribution\n" << endl;
        while (simple.correctNonOrthogonal())
        {
            solve
            (
                fvm::ddt(alpha)
              ==
                n*kprime*Foam::pow(alpha, (n - 1)/n)*Foam::pow(1 - alpha, (n + 1)/n)
            );
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

