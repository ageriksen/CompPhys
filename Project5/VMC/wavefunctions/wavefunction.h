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

        virtual double calculate( const arma::Mat<double> & positions);
        virtual double localEnergy( const arma::Mat<double> & positions);

    protected:
        int m_NParticles, m_NDimensions;
};

inline double Wavefunction::calculate( const arma::Mat<double> & positions)
{ // non-interacting wavefunction
    exit(0);
    return 0.0;
}

inline double Wavefunction::localEnergy( const arma::Mat<double> & positions)
{ // non-interacting wavefunction
    exit(0);
    return 0.0;
}

#endif // end of wavefunction base class
