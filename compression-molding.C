// This file is based on OpenFOAM(R). Simple modification of laplacianFoam.
//
// DISCLAIMER
// This offering is not approved or endorsed by OpenCFD Limited, producer and distributor
// of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM(R) and OpenCFD(R)
// trade marks.
//
// ACKNOWLEDGEMENT
// OPENFOAM®  is a registered trade mark of OpenCFD Limited, producer and distributor of the
// OpenFOAM software via www.openfoam.com.

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
        tsum += runTime.deltaT()/ti;

        volScalarField xi(pos(1 - tsum));
        forAll(xi, i)
        {
            // tsum[i] > 1.0 for the first time
            if (xi.v()[i] >= 1.0 and tinduction.v()[i] < 0.0)
            {
                tinduction[i] = runTime.timeOutputValue();
            }
        }

        volScalarField k(k0*Foam::exp(-Ea/R/T));

        Info<< "\nCalculating cure rate distribution\n" << endl;

        // If tsum in the cell is less the 1.0, alpha == 0, otherwise
        // it is calculated using formula
        alpha =
            xi*k*Foam::pow(runTime - tinduction, n)/(1 + k*Foam::pow(runTime - tinduction, n));

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

