#include "storage.h"
#include "vmcsystem.h"
#include "wavefunctions/wavefunction.h"
#include "wavefunctions/trialwf1naive.h"
#include "wavefunctions/trialwf1full.h"
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
    // System variables
    //------------------------------------------------------
    int NParticles = 2;
    int NDimensions = 3;
    // VMC variables
    int MCCycles = 1e6;
    //------------------------------------------------------
    // Imported Variables
    //------------------------------------------------------
    string baseFileName = "../data/testrun";
    vector<double> parameters(2); // setting up parameter vector to have apropriate size
    vector<string> parameterNames(2);
    vector<double> omegaVec(3);
    vector<double> alphaVec(3);
    vector<double> betaVec(3);
    Wavefunction *WF;
    //

    //if 2nd argument provided, it's the name of the file to store values
    if( numberOfArguments > 2 )
    {
        cout << "test filename triggered" << endl;
        baseFileName = argumentList[2];
        cout << "base filename: " << baseFileName << endl;
    }

    //if 3rd argument provided, choose which wavefunction to use and which loops to run
    if( numberOfArguments > 3 )
    {
        cout << "test wavefunction triggered " << endl;
        int value = atoi(argumentList[3]);
        if( value == 0 )
        {
            cout << "naive hamiltonian" << endl;
            WF = new trialWF1Naive( NParticles, NDimensions );
        }
        else if( value == 1 )
        {
            cout << "full hamiltonian." << endl;
            WF = new trialWF1Full( NParticles, NDimensions );
        }
    }

    //if 4thd argument provided, it's the name of the file with the omega array
    if( numberOfArguments > 4 )
    {
        cout << "test omega triggered" << endl;
        omegaVec.resize(0);
        parameters.resize(0);
        std::string omegaFilename = argumentList[4];
        storage omegaStore(omegaFilename);
        omegaStore.in(&omegaVec);
        parameters.push_back(0.0);
        parameterNames.push_back("omega");
        cout << "size of omegaVec is: " << omegaVec.size() << endl;
        omegaStore.oClose();
        for( auto i: omegaVec )
        {
            cout << i << endl;
        }
    }

    // if 5th argument provided, set alpha start, stop and step from provided file
    if( numberOfArguments > 5 )
    {
        cout << "test alpha triggered" << endl;
        alphaVec.resize(0);
        std::string alphaFilename = argumentList[5];
        storage alphaStore(alphaFilename);
        alphaStore.in(&alphaVec);
        parameters.push_back(0.0);
        parameterNames.push_back("alpha");
        cout << "size of alphaVec is: " << alphaVec.size() << endl;
        alphaStore.oClose();
        for( auto i: alphaVec )
        {
            cout << i << endl;
        }
    }

    // if 5th argument provided, set beta start, stop and step from provided file
    if( numberOfArguments > 6 )
    {
        cout << "test beta triggered" << endl;
        betaVec.resize(0);
        std::string betaFilename = argumentList[6];
        storage betaStore(betaFilename);
        betaStore.in(&betaVec);
        parameters.push_back(0.0);
        parameterNames.push_back("beta");
    }

    string headLine;
    for( unsigned int i = 1; i < parameterNames.size(); i++ )
    {
        cout << "parameterName element: " << parameterNames[i] << endl;
        headLine += parameterNames[i]+" ";
    }
    headLine += "E Variance Ratio";
    cout << "headline set to " << headLine << endl;

    //------------------------------------------------------
    // Program timer
    //------------------------------------------------------
    steady_clock::time_point programStart = steady_clock::now();



    //
    VMCSystem VMC( NParticles, NDimensions, processors, processRank );
    //------------------------------------------------------
    //WAVEFUNCTION
    //------------------------------------------------------
    //
    VMC.setWaveFunction( WF );


    //------------------------------------------------------
    //  Setting stepLength
    //------------------------------------------------------
    vector<double> stepLength(omegaVec.size());
    for( unsigned int index= 0; index < omegaVec.size(); index++ )
    {
        std::cout
            << "\n------------------------------------------------------\n"
            << "\n omega = " << omegaVec[index] << "\n";
        parameters[0] = omegaVec[index];
        parameters[1] = 1;
        stepLength[index] = VMC.stepFinder( parameters );
    }

    //------------------------------------------------------
    // looping over variables
    //------------------------------------------------------
    cout << headLine << endl;
    string fileName;
    for( unsigned int index= 0; index < omegaVec.size(); index++ )
    {
        std::cout
            << "\n------------------------------------------------------\n";
        fileName = baseFileName + "Omega" + std::to_string(omegaVec[index]) + ".dat";

        storage storageFile(fileName);
        storageFile.out();
        storageFile.dat(headLine);

        for( double alpha = alphaVec[0]; alpha < alphaVec[1]; alpha += alphaVec[2])
        {
            //running variationaVecl MC
            parameters[0] = omegaVec[index];
            parameters[1] = alpha;
            WF -> setParameters( parameters );
            VMC.runVMC( MCCycles, stepLength[index] );

            // storing run
            if( processRank == 0 )
            {
                storageFile.lineAdd( std::to_string(alpha) );
                storageFile.lineAdd( std::to_string(VMC.energy()/double(MCCycles)) );
                storageFile.lineAdd( std::to_string( ( VMC.energySquared() - (VMC.energy()*VMC.energy()/MCCycles) )/double(MCCycles) ) );
                storageFile.lineAdd( std::to_string( VMC.ratio() ) );
                cout << storageFile.getLine() << endl;
                storageFile.dat();
                storageFile.lineClean();
            }
        } // end for alpha
        storageFile.oClose();
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
        //storageFile.oClose();
    }

    //------------------------------------------------------
    // MPI namespace cancel:
    //------------------------------------------------------
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();

    return 0;
} // end main
