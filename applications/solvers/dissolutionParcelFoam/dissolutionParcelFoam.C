/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | User-space dissolution solver for OpenFOAM 9
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "kinematicMomentumTransportModel.H"
#include "pisoControl.H"
#include "pressureReference.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "DissolutionCloud.H"

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"
        #include "setDeltaT.H"

        MRF.update();

        {
            fvModels.correct();

            #include "UEqn.H"

            while (piso.correct())
            {
                #include "pEqn.H"
            }
        }

        dissolutionCloud.evolve();

        #include "CEqn.H"

        laminarTransport.correct();
        turbulence->correct();

        Info<< "Active parcels          = " << dissolutionCloud.size() << nl
            << "Total parcel mass [kg]  = "
            << dissolutionCloud.totalParticleMass() << nl
            << "Integrated scalar mass  = "
            << fvc::domainIntegrate(C).value() << nl << endl;

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
