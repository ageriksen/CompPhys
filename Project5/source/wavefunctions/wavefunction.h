#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <armadillo>
#include <cmath>
#include <string>

class Wavefunction
{
    public:
        Wavefunction(int NParticles, int NDimensions):
            m_NParticles(NParticles), m_NDimensions(NDimensions)
        {}

        virtual void setParameters( arma::Col<double> parameters ) = 0;
        virtual double powers( const arma::Mat<double> & positions) = 0;
        virtual double localEnergy( const arma::Mat<double> & positions) = 0;
        //virtual std::string save() = 0;

    protected:
        int m_NParticles, m_NDimensions;
};

#endif // end of wavefunction base class
