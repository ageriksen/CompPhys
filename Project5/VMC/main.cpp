#include "store.h"
#include "vmcsystem.h"
#include "wavefunctions/twoparticlenoninteractingwf.h"
#include <mpi.h>
#include <iostream>

using namespace std;
using chrono::steady_clock;
using chrono::duration_cast;
using chrono::duration;

int main( int numberOfArguments, char *cmdLineArguments[])
{

    //MPI initialization
    int processors, processRank;
    MPI_Init( &numberOfArguments, &cmdLineArguments );
    MPI_Comm_size( MPI_COMM_WORLD, &processors );
    MPI_Comm_rank( MPI_COMM_WORLD, &processRank );

    // Program timer
    steady_clock::time_point programStart = steady_clock::now();

    // System variables
    int NParticles = 2;
    int NDimensions = 3;
    double omega = 1.0;
    // VMC variables
    int MCCycles = 1e5;
    double steplength = 1.0; // should give ~50% acceptance
    //
    VMCSystem VMC( NParticles, NDimensions, processors, processRank );
    TwoParticleNonInteractingWF WF( NParticles, NDimensions );
    VMC.setWaveFunction( &WF );

    // Store object
    string fileName;
    cout << "please provide a filename for storage" << endl;
    cin >> fileName;
    store alphaFile( fileName );
    alphaFile.open();

    string headline = "alpha E Variance {Acceptance ratio}";
    alphaFile.dat(headline);


    // Varying wavefunction parameters
    double alphaMin = 0.8;
    double alphaMax = 1.2;
    double alphaStep = 0.05;
    for( double alpha = alphaMin; alpha < alphaMax; alpha += alphaStep )
    {
        if( processRank == 0 )
        {
            cout << "\n alpha = " << alpha << "\n";
        }
        //running variational MC
        WF.setParameters( omega, alpha );
        VMC.runVMC( MCCycles, steplength );

        // storing run
        if( processRank == 0 )
        {
            alphaFile.lineAdd( to_string(alpha) );
            alphaFile.lineAdd( to_string(VMC.energy()/double(MCCycles)) );
            alphaFile.lineAdd( to_string( ( VMC.energySquared() - VMC.energy() )/double(MCCycles) ) );
            alphaFile.lineAdd( to_string( VMC.ratio() ) );
            alphaFile.dat();
            alphaFile.lineClean();
        }
    }

    // printing runtime
    duration<double> programTime = duration_cast<duration<double>>(steady_clock::now() - programStart);

    if( processRank == 0 )
    {
        cout << "runtime complete. time used: \n"
             << double(programTime.count())/3600.0 << " hrs ("
             << programTime.count() << ")" << endl;
        alphaFile.close();
    }

    // MPI namespace cancel:
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();

    return 0;
} // end main
