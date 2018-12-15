#ifndef VMCSYSTEM_H
#define VMCSYSTEM_H

#include <vector>
#include "wavefunctions/wavefunction.h"

using std::vector;

class VMCSystem
{
    public:
        VMCSystem( int NParticles, int NDimensions, int processors, int rank):
            m_NParticles(NParticles), m_NDimensions(NDimensions),
            m_processors(processors), m_rank(rank)
        {
            vector<double> dimensions(m_NDimensions, 0.0);
            m_positionsOld.resize(m_NParticles, dimensions);
            m_positionsNew.resize(m_NParticles, dimensions);
        }
        //
        void setWaveFunction( Wavefunction * WF )
        {
            m_WF = WF;
        }
        //
        void runVMC( int MCCycles, double steplength );
        double stepFinder( vector<double> param );
        void clean()
        {
        m_energy = 0;
        m_energySquared = 0;
        m_oldWaveFunction = 0;
        m_newWaveFunction = 0;
        }

        // getters
        double energy() { return m_energy; }
        double energySquared() { return m_energySquared; }
        double ratio() { return m_acceptRatio; }

    private:
        //----------------------------------------
        // System variables
        //----------------------------------------
        int m_NParticles, m_NDimensions;
        vector< vector<double> > m_positionsOld, m_positionsNew;
        // Wavefunction
        Wavefunction * m_WF = nullptr;
        double m_oldWaveFunction = 0;
        double m_newWaveFunction = 0;
        // Parallellization variables
        int m_processors, m_rank;
        // MC variables
        int m_MCCycles;
        double m_acceptRatio;
        // resultant variables
        double m_energy = 0;
        double m_energySquared = 0;
}; // end of system class

#endif // end of VMCSYSTEM_H
