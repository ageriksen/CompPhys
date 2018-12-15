#include "vmcsystem.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>
#include <mpi.h>

//------------------------------------------------------
//  VRUN MC
//------------------------------------------------------
void VMCSystem::runVMC( int MCCycles, double steplength )
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

    //------------------------------------------------------
    // initial wavefunction exponent
    //------------------------------------------------------
    m_oldWaveFunction = m_WF -> powers(m_positionsOld);

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
            m_newWaveFunction = m_WF -> powers(m_positionsNew);

            // updates ratio:
            exponent = 2.0*m_newWaveFunction - 2.0*m_oldWaveFunction;
            ratio = std::exp(exponent);

            if( acceptanceDistribution( engine ) <= ratio )
            {// accepting suggestion
                m_oldWaveFunction = m_newWaveFunction;
                for( int dimension = 0; dimension < m_NDimensions; dimension ++ )
                {
                    m_positionsOld[particle][dimension] = m_positionsNew[particle][dimension];
                }
                // Update acceptance:
                acceptance++;
            }
            else
            {//reset parameters
                m_newWaveFunction = m_oldWaveFunction;
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


    //------------------------------------------------------
    // Gathering data from processors
    //------------------------------------------------------
    double tmpEnergy = 0;
    double tmpEnergySquared = 0;
    double tmpAcceptRatio = 0;

    MPI_Reduce( &m_energy,          &tmpEnergy,         1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce( &m_energySquared,   &tmpEnergySquared,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce( &m_acceptRatio,     &tmpAcceptRatio,    1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


} // end runVMC


//------------------------------------------------------
//  STEP FINDER
//------------------------------------------------------
double VMCSystem::stepFinder( vector<double> param )
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
        this -> m_WF -> setParameters( param );
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
