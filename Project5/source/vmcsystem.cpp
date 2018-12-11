#include "vmcsystem.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <mpi.h>

using std::cout;
using std::endl;
using std::setprecision;

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
    double localEnergy = 0;
    m_energy = 0;
    m_energySquared = 0;
    //reset acceptance as well
    m_acceptRatio = 0;

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
    unsigned long long acceptance = 0;
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

            if( acceptanceDistribution( engine ) <= ratio )
            {// accepting suggestion
                m_oldWaveFunction = m_newWaveFunction;
                for( int dimension = 0; dimension < m_NDimensions; dimension ++ )
                {
                    m_positionsOld(particle, dimension) = m_positionsNew(particle, dimension);
                }
                // Update acceptance:
                acceptance++;
            }
            else
            {//reset parameters
                m_newWaveFunction = m_oldWaveFunction;
                for( int dimension = 0; dimension < m_NDimensions; dimension ++ )
                {
                    m_positionsNew(particle, dimension) = m_positionsOld(particle, dimension);
                }
            }
        } // end of metropolis per particle

        // New local energy
        localEnergy = m_WF -> localEnergy(m_positionsNew);
        m_energy += localEnergy;
        m_energySquared += localEnergy*localEnergy;


    } // end of MCCycles
    m_acceptRatio = double(acceptance)/double(m_NParticles*m_MCCycles);


    // Gathering data from processors
    double tmpEnergy = 0;
    double tmpEnergySquared = 0;
    double tmpAcceptRatio = 0;

    MPI_Reduce( &m_energy,          &tmpEnergy,         1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce( &m_energySquared,   &tmpEnergySquared,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce( &m_acceptRatio,     &tmpAcceptRatio,    1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // printout
    if( m_rank == 0 )
    {
        cout << setprecision(16) << "Mean energy:       " << m_energy/double(MCCycles) << "\n"
                                 << "Variance(Energy):  " << (m_energySquared - (m_energy*m_energy/MCCycles)) / double(MCCycles) << "\n"
                                 << "Acceptance ratio:  " << m_acceptRatio
             << endl;
    }


} // end runVMC
