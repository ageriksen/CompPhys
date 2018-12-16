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
    m_distance = 0;
    m_kinetic = 0;
    m_potential = 0;

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
        m_distance += m_wF->distance();
        m_kinetic += m_WF->kinetic();
        m_potential += m_WF->potential();


    } // end of MCCycles
    m_acceptRatio = double(acceptance)/double(m_NParticles*m_MCCycles);


    //------------------------------------------------------
    // Gathering data from processors
    //------------------------------------------------------
    double tmpEnergy = 0;
    double tmpEnergySquared = 0;
    double tmpAcceptRatio = 0;
    double tmpDistance = 0;
    double tmpKinetic = 0;
    double tmpPotential = 0;

    MPI_Reduce( &m_energy,          &tmpEnergy,         1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce( &m_energySquared,   &tmpEnergySquared,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce( &m_acceptRatio,     &tmpAcceptRatio,    1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce( &m_distance,        &tmpDistance,       1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce( &m_kinetic,         &tmpKinetic,        1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce( &m_potential,       &tmpPotential,      1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    m_energy = tmpEnergy;
    m_energySquared = tmpEnergySquared;
    m_acceptRatio = tmpAcceptRatio/m_processors;
    m_distance = tmpDistance/double(m_MCCyclesFull*m_processors);
    m_kinetic = tmpKinetic/double(m_MCCyclesFull*m_processors);
    m_potential= tmpPotential/double(m_MCCyclesFull*m_processors);

    m_variance = ( m_energySquared - m_energy*m_energy/double(m_MCCyclesFull*m_processors))/double(m_MCCyclesFull*m_processors);
    m_energyMean = m_energy/double(m_MCCyclesFull*m_processors);


} // end runVMC


//------------------------------------------------------
//  STEP FINDER
//------------------------------------------------------
double VMC::stepFinder( vector<double> param )
{
    // variable setup
    double trialStep;
    int quickMC = 1e5; // only need acceptance ratio, not exact results
    double acceptanceMin = 0.48; // lower bound on acceptance ratio
    double acceptanceMax = 0.53; // upper bound on acceptance ratio
    double currentAccept;
    double deltaMin = 0.01/double(param[0]); // min delta value tested and step
    double deltaMax = 10/double(param[0]);  // max delta value tested
    std::cout << "min and step: " << deltaMin << ", max: " << deltaMax <<  std::endl;

    //testing steplengths
    for( double delta = deltaMin; delta < deltaMax; delta += deltaMin )
    {
        m_WF -> setParameters( param );
        this -> runVMC( quickMC, delta );

        //std::cout << "delta, acceptratio: " << delta << ", " << m_acceptRatio << "\n";

        if(  m_acceptRatio > acceptanceMin && m_acceptRatio < acceptanceMax)
        {
            currentAccept = m_acceptRatio;
            trialStep = delta;
            std::cout << "steplength: " << trialStep << "\n"
                      << "acceptRatio: " << currentAccept << std::endl;
            clean();
            return trialStep;
        }
    }
    if( m_rank == 0 )
    {
        std::cout << "could not find a good steplength" << std::endl;
//        exit(0);
    }
    clean();
    return 5;
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
    cout << "entered alpha0" << endl;
    double tmpAlpha;
    double varianceMin;
    double delta;
    string tmpLine;
    vector<double> param(2);
    vector<double> alphaMin;
    for( unsigned int omega = 0; omega < parameters[0].size(); omega++ )
    {
        delta = parameters[1][omega];
        varianceMin = 10; // pretty sure the variance should dip WAY below this.
        vector<string> line;
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
                tmpLine = to_string(parameters[2][alpha]) + " "
                        + to_string(m_energyMean) + " "
                        + to_string(m_variance) + " "
                        + to_string(m_acceptRatio);
                line.push_back(tmpLine);
                //cout << tmpLine << endl;
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
            variableSaver.out();
            cout << "-------------------------------------" << endl;
            for( string element: line )
            {
                variableSaver.dat(element);
                cout << element << endl;
            }
            variableSaver.close();
        }
    }
    cout << "time to save this stuff." << endl;
    storage alpha0Saver("./resources/alpha0.dat");
    alpha0Saver.out();
    for( double alpha0: alphaMin )
    {
        alpha0Saver.dat(to_string(alpha0));
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
