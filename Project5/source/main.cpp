#include "storage.h"
#include "vmc.h"
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
    vector< vector<double> > parameters; //
    vector<string> parameterNames;
    vector<double> omegaVec;
    vector<double> stepLength;
    vector<double> alphaVec;
    vector<double> betaVec;
    Wavefunction *WF;
    //
    //-----------------------------
    // establishing base file name
    //-----------------------------
    baseFileName = argumentList[2];
    cout << "base filename: " << baseFileName << endl;

    //-----------------------------
    // deciding on wavefunction type
    //-----------------------------
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

    //-----------------------------
    // getting omega values
    //-----------------------------
    string omegaFilename = argumentList[4];
    storage omegaFinder(omegaFilename);
    omegaFinder.in(&omegaVec);
    omegaFinder.close();
    stepLength.resize(omegaVec.size(), 0 );

    parameters.push_back(omegaVec);
    cout << "size of omegaVec is: " << parameters[0].size() << endl;
    for( double i: parameters[0] )
    {
        cout << i << endl;
    }

    //-----------------------------
    // getting alpha values
    //-----------------------------
    std::string alphaFilename = argumentList[5];
    storage alphaFinder(alphaFilename);
    alphaFinder.in(&alphaVec);
    alphaFinder.close();

    parameters.push_back(alphaVec);
    parameterNames.push_back("alpha");
    if( processRank == 0 )
    {
        cout << "size of alphaVec is: " << alphaVec.size() << endl;
        for( double element : alphaVec )
        {
            cout << element << endl;
        }
        cout << "size of parameters is: " << parameters[2].size() << endl;
        for( double i: parameters[2])
        {
            cout << i << endl;
        }
    }

    //-----------------------------
    // if 5th argument provided, set
    // beta start, stop and step from
    // provided file
    //-----------------------------
    if( numberOfArguments > 6 )
    {
        cout << "test beta triggered" << endl;
        std::string betaFilename = argumentList[6];
        storage betaStore(betaFilename);
        betaStore.in(&betaVec);
        betaStore.close();

        parameters.push_back(betaVec);
        parameterNames.push_back("beta");
        cout << "size of betaVec is: " << betaVec.size() << endl;
        for( double i: betaVec )
        {
            cout << i << endl;
        }
    }

    //------------------------------------------------------
    // Program timer
    //------------------------------------------------------
    steady_clock::time_point programStart = steady_clock::now();



    //------------------------------------------------------
    // declaring VMC
    //------------------------------------------------------
    VMC vmc( NParticles, NDimensions, processors, processRank );
    //------------------------------------------------------
    //WAVEFUNCTION
    //------------------------------------------------------
    //
    vmc.setWaveFunction( WF );


    //------------------------------------------------------
    //  Setting stepLength
    //------------------------------------------------------
    vector<double> param(2);
    for( unsigned int index= 0; index < parameters[0].size(); index++ )
    {
        cout
            << "\n------------------------------------------------------\n"
            << "\n omega = " << parameters[0][index] << "\n";
        param[0] = parameters[0][index];
        param[1] = 1;
        cout << "entering stepFinder" << endl;
        parameters[1][index] = vmc.stepFinder( param );
    }

    //------------------------------------------------------
    // Finding ideal alpha, storing variables
    //------------------------------------------------------
    vector<double> alpha0 = vmc.alpha0(parameters, MCCycles, baseFileName );

    //------------------------------------------------------
    // printing runtime
    //------------------------------------------------------
    duration<double> programTime = duration_cast<duration<double>>(steady_clock::now() - programStart);

    if( processRank == 0 )
    {
        std::cout << "runtime complete. time used: \n"
             << double(programTime.count())/3600.0 << " hrs ("
             << programTime.count() << ")" << std::endl;
        //storageFile.close();
    }

    //------------------------------------------------------
    // MPI namespace cancel:
    //------------------------------------------------------
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();

    return 0;
} // end main
