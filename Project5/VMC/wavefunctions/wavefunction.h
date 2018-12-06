#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <armadillo>

class Wavefunction
{
    public:
        Wavefunction(int NParticles, int NDimensions):
            m_NParticles(NParticles), m_NDimensions(NDimensions)
    {}
        virtual ~Wavefunction() {}

        virtual double ratio( const arma::Mat<double> & Position);
        virtual double localEnergy( const arma::Mat<double> & Position);

    protected:
        int m_NParticles, m_NDimensions;
};

inline double Wavefunction::ratio( const arma::Mat<double> & Position)
{ // non-interacting wavefunction
    exit(0);
    return 0.0;
}

inline double Wavefunction::localEnergy( const arma::Mat<double> & Position)
{ // non-interacting wavefunction
    exit(0);
    return 0.0;
}

#endif // end of wavefunction base class
