#include "store.h"
#include "vmcsystem.h"
#include "wavefunctions/wavefunction.h"
#include "wavefunctions/trialwf1naive.h"
#include <mpi.h>
#include <iostream>

//------------------------------------------------------
//  NAMESPACES
//------------------------------------------------------
using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::duration;

//------------------------------------------------------
//  MAIN
//------------------------------------------------------
int main( int numberOfArguments, char *cmdLineArguments[])
{

    //------------------------------------------------------
    //MPI initialization
    //------------------------------------------------------
    int processors, processRank;
    MPI_Init( &numberOfArguments, &cmdLineArguments );
    MPI_Comm_size( MPI_COMM_WORLD, &processors );
    MPI_Comm_rank( MPI_COMM_WORLD, &processRank );

    //------------------------------------------------------
    // Program timer
    //------------------------------------------------------
    steady_clock::time_point programStart = steady_clock::now();

    //------------------------------------------------------
    // System variables
    //------------------------------------------------------
    int NParticles = 2;
    int NDimensions = 3;
    double omega = 1.0;
    // VMC variables
    int MCCycles = 1e6;
    double stepLength = 1.0; // should give ~50% acceptance
    //
    VMCSystem VMC( NParticles, NDimensions, processors, processRank );
    //------------------------------------------------------
    //WAVEFUNCTION
    //------------------------------------------------------
    Wavefunction *WF;
    //
    WF = new trialWF1Naive( NParticles, NDimensions );
    //
    VMC.setWaveFunction( WF );

    //------------------------------------------------------
    // Storage object
    //------------------------------------------------------
    string fileName;
    std::cout << "please provide a filename for storage" << std::endl;
    std::cin >> fileName;
    store alphaFile( fileName );
    alphaFile.open();
    //
    string headline = "alpha E Variance {Acceptance ratio}";
    alphaFile.dat(headline);

    //------------------------------------------------------
    //  RUNNINGN VMC
    //------------------------------------------------------
    // Varying wavefunction parameters
    double alphaMin = 0.8;
    double alphaMax = 1.2;
    double alphaStep = 0.05;
    for( double alpha = alphaMin; alpha < alphaMax; alpha += alphaStep )
    {
        if( processRank == 0 )
        {
            std::cout << "\n alpha = " << alpha << "\n";
        }
        //running variational MC
        arma::Col<double> parameters = ( parameters );
        WF -> setParameters( parameters );
        VMC.runVMC( MCCycles, stepLength );

        // storing run
        if( processRank == 0 )
        {
            alphaFile.lineAdd( std::to_string(alpha) );
            alphaFile.lineAdd( std::to_string(VMC.energy()/double(MCCycles)) );
            alphaFile.lineAdd( std::to_string( ( VMC.energySquared() - (VMC.energy()*VMC.energy()/MCCycles) )/double(MCCycles) ) );
            alphaFile.lineAdd( std::to_string( VMC.ratio() ) );
            alphaFile.dat();
            alphaFile.lineClean();
        }
    }


    //------------------------------------------------------
    // printing runtime
    //------------------------------------------------------
    duration<double> programTime = duration_cast<duration<double>>(steady_clock::now() - programStart);

    if( processRank == 0 )
    {
        std::cout << "runtime complete. time used: \n"
             << double(programTime.count())/3600.0 << " hrs ("
             << programTime.count() << ")" << std::endl;
        alphaFile.close();
    }

    //------------------------------------------------------
    // MPI namespace cancel:
    //------------------------------------------------------
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();

    return 0;
} // end main
