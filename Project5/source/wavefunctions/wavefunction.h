#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <armadillo>
#include <cmath>

class Wavefunction
{
    public:
        Wavefunction(int NParticles, int NDimensions):
            m_NParticles(NParticles), m_NDimensions(NDimensions)
        {}
        virtual ~Wavefunction() {}

        virtual double powers( const arma::Mat<double> & positions) = 0;
        virtual double localEnergy( const arma::Mat<double> & positions) = 0;

    protected:
        int m_NParticles, m_NDimensions;
};

#endif // end of wavefunction base class
