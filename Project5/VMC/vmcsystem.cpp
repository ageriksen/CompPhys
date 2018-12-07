#include "vmcsystem.h"

void VMCSystem::runVMC( int MCCycles, double steplength )
{
    /*
     * Implementation of the variational Monte Carlo method for finding the
     * optimal trial wavefunction of a given family.
     */

    //initialize ratio variable
    double exponent = 0;
    double ratio = 0;

    //ensure energies to be stored are nulled out before the start.
    m_energy = 0;
    m_energySquared = 0;
    //reset acceptance as well
    m_acceptanceCounter = 0;

    // RNG setup
    std::random_device rd;
    std::mt19937_64 engine(rd() + m_rank);
    std::uniform_real_distribution<double> uniformDistribution(-1, 1);
    std::uniform_real_distribution<double> acceptanceDistribution(0, 1);

    // initial value for system
    for( int particle = 0; particle < m_NParticles; particle ++ )
    {
        for( int dimension = 0; dimension < m_NDimensions; dimension ++ )
        {
            m_positionsOld(particle, dimension) = uniformDistribution( engine );
            m_positionsNew(particle, dimension) = m_positionsOld(particle, dimension);
        }
    }

    // MCCycles split amongst processors:
    m_MCCycles = MCCycles / m_processors;

    // initial wavefunction exponent
    m_oldWaveFunction = m_WF -> powers(m_positionsOld);

    //----------------------------------------
    //  MC simulation over total cycles
    //----------------------------------------
    for( int cycle = 0; cycle < m_MCCycles; cycle ++ )
    {

        // runs each particle through the metropolis sampling rule
        for( int particle = 0; particle < m_NParticles; particle ++ )
        {

            //Suggest update
            for( int dimension = 0; dimension < m_NDimensions; dimension ++ )
            {
                m_positionsNew(particle, dimension) = m_positionsOld(particle, dimension)
                                                    + steplength*uniformDistribution( engine );
            } // end update suggestion
            // and finds new wavefunction exponent
            m_newWaveFunction = m_WF -> powers(m_positionsNew);

            // updates ratio:
            exponent = 2.0*m_newWaveFunction - 2.0*m_oldWaveFunction;
            ratio = std::exp(exponent);

        } // end of metropolis per particle
    } // end of MCCycles


} // end runVMC
