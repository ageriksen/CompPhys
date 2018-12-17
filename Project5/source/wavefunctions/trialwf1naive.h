#ifndef TRIALWF1NAIVE_H
#define TRIALWF1NAIVE_H

#include "wavefunction.h"
#include <iostream>

class trialWF1Naive: public Wavefunction
{
//------------------------------------------------------
    public:
        trialWF1Naive( int NParticles, int NDimensions );

        void setParameters( vector<double> parameters )
        {
            m_omega = parameters[0];
            m_alpha = parameters[1];
        }

        double sumSquares( const vector< vector<double> > & positions );

        double powers( const vector< vector<double> > & positions );
        double localEnergy( const vector< vector<double> > & positions );

        //getters and setters
        double omega()  { return m_omega; }
        double alpha()  { return m_alpha; }
        double m_distance() { return 0; }
        double m_kinetic() { return 0; }
        double m_potential() { return 0; }

//------------------------------------------------------
    private:
        double m_omega = 0;
        double m_alpha = 0;
}; //

//------------------------------------------------------
inline double trialWF1Naive::powers
(
 const vector< vector<double> > & positions
)
{ // powers of state wavefunction.
    return -0.5*m_alpha*m_omega*sumSquares(positions);
}

//------------------------------------------------------
inline double trialWF1Naive::localEnergy
(
 const vector< vector<double> > & positions
)
{ // local energy of naive trial wavefunction
    return 0.5*m_omega*m_omega
        *sumSquares(positions)
        *( 1 - m_alpha*m_alpha )
        + 3.0*m_alpha*m_omega;
}

#endif // TRIALWF1NAIVE_H
