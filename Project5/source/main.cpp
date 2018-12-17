#include "storage.h"
#include "vmc.h"
#include "wavefunctions/wavefunction.h"
#include "wavefunctions/trialwf1naive.h"
#include "wavefunctions/trialwf1full.h"
#include "wavefunctions/trialwf2.h"
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
using std::to_string;

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
    double r12 = 0;
    // VMC variables
    int MCCycles = 1e6;

    //------------------------------------------------------
    // Imported Variables
    //------------------------------------------------------
    string baseFileName = "../data/testrun";
    vector< vector<double> > parameters; //
    vector<string> parameterNames;
    vector<double> omegaVec;
    vector<double> deltaVec;
    vector<double> alphaVec;
    vector<double> betaVec;
    double alpha0;
    Wavefunction *WF = nullptr;
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
    else if( value == 2 )
    {
        cout << "trial wavefunction 2 " << endl;
        WF = new trialWF2( NParticles, NDimensions );
    }
    else
    {
        cout << "Please provide a wavefunction to run for" << endl;
        return 1;
    }

    //------------------------------------------------------
    // declaring VMC
    //------------------------------------------------------
    VMC vmc( NParticles, NDimensions, processors, processRank );
    vmc.setWaveFunction( WF );

    //-----------------------------
    // getting omega values
    //-----------------------------
    string omegaFilename = argumentList[4];
    storage omegaFinder(omegaFilename);
    omegaFinder.in(&omegaVec);
    omegaFinder.close();
    deltaVec.resize(omegaVec.size(), 0 );

    parameters.push_back(omegaVec);
    parameters.push_back(omegaVec); // second one is for the steplength
    if( processRank == 0 )
    {
        cout << "size of omegaVec is: " << parameters[0].size() << endl;
        for( double i: parameters[0] )
        {
            cout << i << endl;
        }
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
        cout << "size of alphaVec is: " << parameters[2].size() << endl;
        for( double i: parameters[2])
        {
            cout << i << endl;
        }
    }
    //------------------------------------------------------
    //  Setting stepLength
    //------------------------------------------------------
    // if test accepts, read from file, else find stepLength
    if( numberOfArguments > 6 )
    {
        cout << "stepLength test triggered" << endl;
        string deltaFileName = argumentList[6];
        storage deltaFinder(deltaFileName);
        deltaFinder.in(&deltaVec);
        deltaFinder.close();

        cout << "size of deltaVec is: " << deltaVec.size() << endl;
        parameters[1] = deltaVec;
        for( double element: parameters[1] )
        {
            cout << element << endl;
        }
    }
    else
    {
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
        storage deltaSaver("./resources/delta.dat");
        deltaSaver.out();
        for( double delta: parameters[1] )
        {
            deltaSaver.dat( to_string(delta) );
        }
        deltaSaver.close();
    }

    //-----------------------------
    // if 5th argument provided, set
    // beta start, stop and step from
    // provided file
    //-----------------------------
    if( numberOfArguments > 7 )
    {
        cout << "test beta triggered" << endl;
        std::string betaFilename = argumentList[7];
        storage betaFinder(betaFilename);
        betaFinder.in(&betaVec);
        betaFinder.close();

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
    //WAVEFUNCTION
    //------------------------------------------------------
    //
    vmc.setWaveFunction( WF );



    //------------------------------------------------------
    // Finding ideal alpha, storing variables
    //------------------------------------------------------
    //if test returns true, find minima from file, else run.
    if( numberOfArguments > 8 )
    {
        cout << "test alpha0 triggered" << endl;
        string alpha0FileName = argumentList[8];
        vector<double> tmpVec;
        storage alpha0Finder(alpha0FileName);
        alpha0Finder.in(&tmpVec);
        alpha0Finder.close();
        alpha0 = tmpVec[0];
        cout << "alpha0 = " << alpha0 << endl;
    }
    else
    {
        alpha0 = vmc.alpha0(parameters, MCCycles, baseFileName );
    }


    //------------------------------------------------------
    //  OPTIMIZING ALPHA AND BETA
    //------------------------------------------------------
    // only if beta is read in:
    if( numberOfArguments > 7 )
    {
        unsigned int n = 4;
        vector<int> range(n);
        vector<double> minima(2);
        vector<double> optimized(3);
        for( unsigned int i=0; i<n; i++)
        {
            range[i] = i;
        }

        string omegaName = baseFileName + "optimized.dat";
        storage optimalSaver(omegaName);
        string headLine = "omega alpha beta E variance acceptance";
        optimalSaver.out();
        optimalSaver.dat(headLine);

        for( unsigned index = 0; index < parameters[0].size(); index++ )
        {
            cout << "optimizing. omega = " << parameters[0][index] << endl;

            optimized = vmc.optimize( WF, parameters, range, index, alpha0 );

            cout << "optimal values for omega " << optimized[0] << " are alpha=" << optimized[1]
                 << " and beta=" << optimized[2] << endl;
            WF -> setParameters( optimized );
            vmc.runVMC( MCCycles*100, vmc.stepFinder( optimized ), r12 );

            optimalSaver.lineAdd(
                    to_string(optimized[0]) + " "
                  + to_string(optimized[1]) + " "
                  + to_string(optimized[2]) + " "
                  + to_string(vmc.meanEnergy()) + " "
                  + to_string(vmc.variance()) + " "
                  + to_string(r12) + " "
                  + to_string(vmc.ratio())
                                );
            optimalSaver.dat();
        }
        optimalSaver.close();
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
        //storageFile.close();
    }

    //------------------------------------------------------
    // MPI namespace cancel:
    //------------------------------------------------------
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();

    return 0;
} // end main
