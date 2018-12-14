#include "storage.h"
#include "vmcsystem.h"
#include "wavefunctions/wavefunction.h"
#include "wavefunctions/trialwf1naive.h"
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <string>

//------------------------------------------------------
//  NAMESPACES
//------------------------------------------------------
using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::cout;
using std::cin;
using std::endl;
using std::string;
using std::vector;

//------------------------------------------------------
//  MAIN
//------------------------------------------------------
int main( int numberOfArguments, char *argumentList[])
{

    //------------------------------------------------------
    //MPI initialization
    //------------------------------------------------------
    int processors, processRank;
    MPI_Init( &numberOfArguments, &argumentList);
    MPI_Comm_size( MPI_COMM_WORLD, &processors );
    MPI_Comm_rank( MPI_COMM_WORLD, &processRank );

    //------------------------------------------------------
    // Imported Variables
    //------------------------------------------------------
    storage storageFile;
    vector<double> parameters; // setting up parameter vector to have apropriate size
    vector<string> parameterNames;
    vector<double> omegaVec;
    vector<double> alphaVec;
    vector<double> betaVec;
    //if 2nd argument provided, it's the name of the file to store values
    if( numberOfArguments > 2 )
    {
        //string fileName = argumentList[2];
        storageFile.name(argumentList[2]);
    }
    //
    //if 3rd argument provided, it's the name of the file with the omega array
    if( numberOfArguments > 3 )
    {
        std::string omegaFilename = argumentList[3];
        storage omegaStore(omegaFilename);
        omegaStore.in(omegaVec);
        parameters.push_back(0.0);
        parameterNames.push_back("omega");
    }
    // if 4th argument provided, set alpha start, stop and step from provided file
    if( numberOfArguments > 4 )
    {
        std::string alphaFilename = argumentList[4];
        storage alphaStore(alphaFilename);
        alphaStore.in(alphaVec);
        parameters.push_back(0.0);
        parameterNames.push_back("alpha");
    }
    // if 5th argument provided, set beta start, stop and step from provided file
    if( numberOfArguments > 5 )
    {
        std::string betaFilename = argumentList[5];
        storage betaStore(betaFilename);
        betaStore.in(betaVec);
        parameters.push_back(0.0);
        parameterNames.push_back("beta");
    }


    //------------------------------------------------------
    // Program timer
    //------------------------------------------------------
    steady_clock::time_point programStart = steady_clock::now();

    //------------------------------------------------------
    // System variables
    //------------------------------------------------------
    int NParticles = 2;
    int NDimensions = 3;
    //double omega = 1.0;
    // VMC variables
    int MCCycles = 1e6;
    //double stepLength; //= 1.0; // should give ~50% acceptance
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
    //      STORAGE INITIATE
    //------------------------------------------------------



    string headline; //= "alpha E Variance {Acceptance ratio}";
    string headBit;
    int headMembers;
    cout << "number of head members: " << endl;
    cin >> headMembers;
    //cout << "number of elements in headline: " << headMembers << endl;
    for( int count = 0; count < headMembers; count ++ )
    {
        cout << "\nnext head element: ";
        cin >> headBit;
        headline += headBit+" ";
        headBit = "";
    }
    storageFile.dat(headline);
    cout << "\n----------------------------------------------\n"
        << headline << endl;

    //------------------------------------------------------
    //  RUNNING VMC
    //------------------------------------------------------

    //------------------------------------------------------
    //      INPUT VALUES
    //------------------------------------------------------
    int omegaSize;
    cout << "\nomegaSize: \n";
    cin >> omegaSize;
    cout << "omegaSize: " << omegaSize << endl;
    arma::Col<double> omega = arma::zeros(omegaSize);
    arma::Col<double> stepLength = arma::zeros(omega.size());
    for( int index = 0; index < omegaSize; index++ )
    {
        cout << "\nomega("<<index<<"): ";
        cin >>  omega(index);
    }

    for( int index= 0; index < omegaSize; index++ )
    {
        std::cout
            << "\n------------------------------------------------------\n"
            << "\n omega = " << omega(index) << "\n";
        parameters[0] = omega(index);
        parameters[1] = 1;
        stepLength(index) = VMC.stepFinder( parameters );
    }
    std::cout
        << "\n------------------------------------------------------\n";

    cout << headline << endl;
    for( unsigned int index= 0; index < omegaVec.size(); index++ )
    {
        cout << "\n----------------------------------------------------\n";
        for( double alpha = alphaVec[0]; alpha < alphaVec[1]; alpha += alphaVec[2])
        {
//            if( processRank == 0 )
//            {
//                std::cout
//                    << "\n------------------------------------------------------\n"
//                    << " alpha = " << alpha << "\n";
//            }
            //running variational MC
            parameters[0] = omega(index);
            parameters[1] = alpha;
            WF -> setParameters( parameters );
            VMC.runVMC( MCCycles, stepLength(index) );

            // storing run
            if( processRank == 0 )
            {
                storageFile.lineAdd( std::to_string(omega(index)) );
                storageFile.lineAdd( std::to_string(alpha) );
                storageFile.lineAdd( std::to_string(VMC.energy()/double(MCCycles)) );
                storageFile.lineAdd( std::to_string( ( VMC.energySquared() - (VMC.energy()*VMC.energy()/MCCycles) )/double(MCCycles) ) );
                storageFile.lineAdd( std::to_string( VMC.ratio() ) );
                cout << storageFile.getLine() << endl;
                storageFile.dat();
                storageFile.lineClean();
            // printout
               // std::cout << std::setprecision(16)
               //     << "Mean energy:       " << VMC.energy()/double(MCCycles) << "\n"
               //     << "Variance(Energy):  " << (VMC.energySquared() - (VMC.energy()*VMC.energy()/double(MCCycles))) / double(MCCycles) << "\n"
               //     << "Acceptance ratio:  " << VMC.ratio()
               //     << std::endl;
            }

        } // end for alpha
    } // end for omega


    //------------------------------------------------------
    // printing runtime
    //------------------------------------------------------
    duration<double> programTime = duration_cast<duration<double>>(steady_clock::now() - programStart);

    if( processRank == 0 )
    {
        std::cout << "runtime complete. time used: \n"
             << double(programTime.count())/3600.0 << " hrs ("
             << programTime.count() << ")" << std::endl;
        storageFile.oClose();
    }

    //------------------------------------------------------
    // MPI namespace cancel:
    //------------------------------------------------------
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();

    return 0;
} // end main
