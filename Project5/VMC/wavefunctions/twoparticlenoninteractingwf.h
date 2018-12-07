#ifndef TWOPARTICLENONINTERACTINGWF_H
#define TWOPARTICLENONINTERACTINGWF_H

#include "wavefunction.h"

class TwoParticleNonInteractingWF : public Wavefunction
{
    public:
        TwoParticleNonInteractingWF( int NParticles, int NDimensions );

        void setParameters( double omega, double alpha )
        {
            m_omega = omega;
            m_alpha = alpha;
        }

        double sumSquares( const arma::Mat<double> & positions );

        double powers( const arma::Mat<double> & positions );
        double localEnergy( const arma::Mat<double> & positions );

    private:
        double m_omega = 0;
        double m_alpha = 0;
}; //

inline double TwoParticleNonInteractingWF::powers
(
 const arma::Mat<double> & positions
)
{ // powers of state wavefunction.
    return -0.5*m_alpha*m_omega*sumSquares(positions);
}

inline double TwoParticleNonInteractingWF::localEnergy
(
 const arma::Mat<double> & positions
)
{ // local energy of naive trial wavefunction
    return 0.5*m_omega*m_omega
        *sumSquares(positions)
        *( 1 - m_alpha*m_alpha )
        + 3.0*m_alpha*m_omega;
}

#endif // TWOPARTICLENONINTERACTINGWF_H
