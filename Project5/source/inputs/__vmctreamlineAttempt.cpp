#include "vmc.h"
#include "storage.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>
#include <mpi.h>

using std::to_string;
using std::cout;
using std::endl;

//------------------------------------------------------
//  VRUN MC
//------------------------------------------------------
void VMC::runVMC( int MCCycles, double steplength )
{
    /*
     * Implementation of the variational Monte Carlo method for finding the
     * optimal trial wavefunction of a given family.
     */

    //------------------------------------------------------
    //initialize ratio variable
    //------------------------------------------------------
    double exponent = 0;
    double ratio = 0;

    //------------------------------------------------------
    //ensure energies to be stored are nulled out before the start.
    //------------------------------------------------------
    double localEnergy = 0;
    m_energy = 0;
    m_energySquared = 0;
    //reset acceptance as well
    m_acceptRatio = 0;

    //------------------------------------------------------
    // RNG setup
    //------------------------------------------------------
    std::random_device rd;
    std::mt19937_64 engine(rd() + m_rank);
    std::uniform_real_distribution<double> uniformDistribution(-1, 1);
    std::uniform_real_distribution<double> acceptanceDistribution(0, 1);

    //------------------------------------------------------
    // initial value for system
    //------------------------------------------------------
    for( int particle = 0; particle < m_NParticles; particle ++ )
    {
        for( int dimension = 0; dimension < m_NDimensions; dimension ++ )
        {
            m_positionsOld[particle][dimension]
                = uniformDistribution( engine );
            m_positionsNew[particle][dimension]
                = m_positionsOld[particle][dimension];
        }
    }



    //------------------------------------------------------
    // MCCycles split amongst processors:
    //------------------------------------------------------
    m_MCCycles = MCCycles / m_processors;
    m_MCCyclesFull = MCCycles;

    //------------------------------------------------------
    // initial wavefunction exponent
    //------------------------------------------------------
    m_oldExponent = m_WF -> powers(m_positionsOld);

    //------------------------------------------------------
    //  MC simulation over total cycles
    //------------------------------------------------------
    unsigned long long acceptance = 0;
    for( int cycle = 0; cycle < m_MCCycles; cycle ++ )
    {
        // runs each particle through the metropolis sampling rule
        for( int particle = 0; particle < m_NParticles; particle ++ )
        {
            //Suggest update
            for( int dimension = 0; dimension < m_NDimensions; dimension ++ )
            {
                m_positionsNew[particle][dimension] = m_positionsOld[particle][dimension]
                                                    + steplength*uniformDistribution( engine );
            } // end update suggestion
            // and finds new wavefunction exponent
            m_newExponent = m_WF -> powers(m_positionsNew);

            // updates ratio:
            exponent = 2.0*m_newExponent - 2.0*m_oldExponent;
            ratio = std::exp(exponent);

            if( acceptanceDistribution( engine ) <= ratio )
            {// accepting suggestion
                m_oldExponent = m_newExponent;
                for( int dimension = 0; dimension < m_NDimensions; dimension ++ )
                {
                    m_positionsOld[particle][dimension] = m_positionsNew[particle][dimension];
                }
                // Update acceptance:
                acceptance++;
            }
            else
            {//reset parameters
                m_newExponent = m_oldExponent;
                for( int dimension = 0; dimension < m_NDimensions; dimension ++ )
                {
                    m_positionsNew[particle][dimension] = m_positionsOld[particle][dimension];
                }
            }
        } // end of metropolis per particle

        //------------------------------------------------------
        // New local energy
        //------------------------------------------------------
        localEnergy = m_WF -> localEnergy(m_positionsNew);
        m_energy += localEnergy;
        m_energySquared += localEnergy*localEnergy;


    } // end of MCCycles
    m_acceptRatio = double(acceptance)/double(m_NParticles*m_MCCycles);

    if( m_rank == 0 )
    {
        cout << m_acceptRatio << endl;
    }

    //------------------------------------------------------
    // Gathering data from processors
    //------------------------------------------------------
    double tmpEnergy = 0;
    double tmpEnergySquared = 0;
    double tmpAcceptRatio = 0;

    MPI_Reduce( &m_energy,          &tmpEnergy,         1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce( &m_energySquared,   &tmpEnergySquared,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce( &m_acceptRatio,     &tmpAcceptRatio,    1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    m_variance = ( m_energySquared - m_energy*m_energy/double(m_MCCyclesFull))/double(m_MCCyclesFull);
    m_energyMean = m_energy/double(m_MCCyclesFull);


} // end runVMC


//------------------------------------------------------
//  STEP FINDER
//------------------------------------------------------
vector<double> VMC::stepFinder( vector<double> parameters )
{
    // variable setup
    double trialStep;
    int quickMC = 1e5; // only need acceptance ratio, not exact results
    double acceptanceMin = 0.48; // lower bound on acceptance ratio
    double acceptanceMax = 0.53; // upper bound on acceptance ratio
    double currentAccept;

    vector<double> deltaVec;
    vector<double> param(2);
    param[1] = 1;
    storage deltaSaver("./resources/delta.dat");

    Wavefunction *WF;
    WF = new trialWF1Naive( 2, 2 );

    for( double omega: parameters )
    {
        param[0] = omega;
        cout << "omega: " << param[0] << endl;
        double deltaMin = 0.01/double(omega); // min delta value tested and step
        double deltaMax = 10/double(omega);  // max delta value tested
        cout << "min and step: " << deltaMin << ", max: " << deltaMax <<  endl;
        //testing steplengths
        for( double delta = deltaMin; delta < deltaMax; delta += deltaMin )
        {
            WF -> setParameters( param );
            this -> runVMC( quickMC, delta );

            if(  m_acceptRatio > acceptanceMin && m_acceptRatio < acceptanceMax)
            {
                currentAccept = m_acceptRatio;
                trialStep = delta;
                std::cout << "steplength: " << trialStep << "\n"
                          << "acceptRatio: " << currentAccept << std::endl;
                deltaVec.push_back(trialStep);
            }
        }
    }
    for( double delta : deltaVec )
    {
        deltaSaver.dat(to_string(delta));
    }
    deltaSaver.close();
    return deltaVec;
}// end stepFinder

//------------------------------------------------------
// OPTIMIZE
//------------------------------------------------------
vector<double> VMC::optimize
(
    vector< vector<double> > parameters,
    vector<int> range,
    double delta,
    double alpha0
)
{
    m_stepLength = delta;
    m_parameters = parameters;
    vector<double> minima(2);
    minima[0] = alpha0;
    minima[1] = findBeta(alpha0);
    for( unsigned int i = 1; i < range.size(); i++ )
    {
        int alpha = findAlpha( minima[1] );
        int beta = findBeta( minima[0] );
        if( alpha < minima[0] && beta < minima[1] )
        {
            minima[0] = alpha;
            minima[1] = beta;
        }
        else
        {
            break;
        }
    }
    return minima;
} // end optimize

//------------------------------------------------------
// FINDERS
//------------------------------------------------------
vector<double> VMC::alpha0
(
    vector< vector<double> > parameters,
    int MCCycles,
    string baseName
)
{
    double tmpAlpha;
    double varianceMin;
    double delta;
    string tmpLine;
    vector<double> param(2);
    vector<double> alphaMin;
    vector<string> line;
    for( unsigned int omega = 0; omega < parameters[0].size(); omega++ )
    {
        delta = parameters[1][omega];
        varianceMin = 10; // pretty sure the variance should dip WAY below this.
        tmpAlpha = 0;
        tmpLine = "";
        for( unsigned int alpha = 0; alpha < parameters[2].size(); alpha++ )
        {
            param[0] = parameters[0][omega];
            param[1] = parameters[2][alpha];
            m_WF -> setParameters( param );
            runVMC( MCCycles, delta );
            if( m_rank == 0 )
            {
                tmpLine += to_string(parameters[2][alpha]) + " "
                        + to_string(m_energyMean) + " "
                        + to_string(m_variance) + " "
                        + to_string(m_acceptRatio);
                cout << tmpLine << endl;
                line.push_back(tmpLine);
                if( m_variance < varianceMin )
                {
                    tmpAlpha = parameters[2][alpha];
                }
            }

        }

        if( m_rank == 0 )
        {
            alphaMin.push_back(tmpAlpha);
            storage variableSaver(baseName + to_string(omega) + ".dat");
            cout << "-------------------------------------" << endl;
            for( string element: line )
            {
                //cout << element << endl;
                variableSaver.dat(element);
            }
            variableSaver.close();
        }
    }
    return alphaMin;
} // end alpha0

double VMC::findAlpha( double beta )
{
    double varianceMin;
    double alphaMin;
    vector<double> param(2);
    param[1] = beta;
    for( double alpha: m_parameters[1] )
    {
        param[0] = alpha;
        m_WF -> setParameters( param );
        runVMC( m_MCCyclesFull, m_stepLength );
        if( m_variance < varianceMin )
        {
            varianceMin = m_variance;
            alphaMin = param[0];
        }
    }
    return alphaMin;
} // end findAlpha


double VMC::findBeta( double alpha )
{
    double varianceMin;
    double betaMin;
    vector<double> param(2);
    param[0] = alpha;
    for( double beta: m_parameters[2] )
    {
        param[1] = beta;
        m_WF -> setParameters( param );
        runVMC( m_MCCyclesFull, m_stepLength );
        if( m_variance < varianceMin )
        {
            varianceMin = m_variance;
            betaMin = param[1];
        }
    }
    return betaMin;
} // end findBeta
