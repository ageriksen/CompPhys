#ifndef VMC_H
#define VMC_H

#include "wavefunctions/wavefunction.h"
#include <vector>
#include <string>

using std::vector;
using std::string;

class VMC
{
    public:
        VMC( int NParticles, int NDimensions, int processors, int rank):
            m_NParticles(NParticles), m_NDimensions(NDimensions),
            m_processors(processors), m_rank(rank)
        {
            vector<double> dimensions(m_NDimensions, 0.0);
            m_positionsOld.resize(m_NParticles, dimensions);
            m_positionsNew.resize(m_NParticles, dimensions);
        }
        //

        //------------------------------
        // class objects
        //------------------------------
        void runVMC( int MCCycles, double steplength );

        double stepFinder( vector<double> param ); // find optimal step for given param

        vector<double> optimize
        (
            vector< vector<double> > parameters,
            vector<int> range,
            double delta, // stepLength
            double alpha0 // best alpha for trial1
        );

        vector<double> alpha0
        ( //find optimal alpha, trial wavefunction 1
         vector< vector<double> > parameters,
         int MCCycles,
         string baseName
        );
        double findAlpha( double beta  ); // find best alpha given beta
        double findBeta ( double alpha ); // find best beta given alpha

        void clean()
        {// set all private variables to 0
        m_energy = 0;
        m_energySquared = 0;
        m_oldExponent = 0;
        m_newExponent = 0;
        m_stepLength = 0;
        m_variance = 0;
        }

        //------------------------------
        // getters, setters
        //------------------------------
        double energy() { return m_energy; }
        double energySquared() { return m_energySquared; }
        double ratio() { return m_acceptRatio; }
        //
        void setWaveFunction( Wavefunction * WF )
        {
            m_WF = WF;
        }

    private:
        //----------------------------------------
        // System variables
        //----------------------------------------
        int m_NParticles, m_NDimensions;
        vector< vector<double> > m_positionsOld;
        vector< vector<double> > m_positionsNew;
        vector< vector<double> > m_parameters;
        // Wavefunction
        Wavefunction * m_WF = nullptr;
        double m_oldExponent = 0;
        double m_newExponent = 0;
        // Parallellization variables
        int m_processors, m_rank;
        // MC variables
        int m_MCCycles;
        int m_MCCyclesFull;
        double m_acceptRatio;
        // resultant variables
        double m_energy = 0;
        double m_energyMean = 0;
        double m_energySquared = 0;
        double m_variance = 0;
        double m_stepLength = 0;
}; // end VMC class

#endif // end of VMCSYSTEM_H
