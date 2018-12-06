#ifndef TWOPARTICLENONINTERACTINGWF_H
#define TWOPARTICLENONINTERACTINGWF_H

#include "wavefunction.h"

class TwoParticleNonInteractingWF : public Wavefunction
{
    public:
        TwoParticleNonInteractingWF( int NParticles, int NDimensions );

        void setParameters( double omega, double alpha );

        double ratio( const arma::mat<double> & Position );

        double localEnergy( const arma::mat<double> & Position );

    private:
        double m_omega = 0;
        double m_alpha = 0;
}; //
#endif // TWOPARTICLENONINTERACTINGWF_H
