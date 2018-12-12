#ifndef TRIALWF1NAIVE_H
#define TRIALWF1NAIVE_H

#include "wavefunction.h"

class trialWF1Naive: public Wavefunction
{
//------------------------------------------------------
    public:
        trialWF1Naive( int NParticles, int NDimensions );

        void setParameters( arma::Col<double> parameters )
        {
            m_omega = parameters(0);
            m_alpha = parameters(0);
        }

        double sumSquares( const arma::Mat<double> & positions );

        double powers( const arma::Mat<double> & positions );
        double localEnergy( const arma::Mat<double> & positions );

        //getters and setters
        double omega()  { return m_omega; }
        double alpha()  { return m_alpha; }

//------------------------------------------------------
    private:
        double m_alpha = 0;
        double m_omega = 0;
}; //

//------------------------------------------------------
inline double trialWF1Naive::powers
(
 const arma::Mat<double> & positions
)
{ // powers of state wavefunction.
    return -0.5*m_alpha*m_omega*sumSquares(positions);
}

//------------------------------------------------------
inline double trialWF1Naive::localEnergy
(
 const arma::Mat<double> & positions
)
{ // local energy of naive trial wavefunction
    return 0.5*m_omega*m_omega
        *sumSquares(positions)
        *( 1 - m_alpha*m_alpha )
        + 3.0*m_alpha*m_omega;
}

#endif // TRIALWF1NAIVE_H
