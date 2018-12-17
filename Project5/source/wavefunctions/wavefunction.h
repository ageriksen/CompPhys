#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

//#include <armadillo>
#include <vector>
#include <cmath>
#include <string>

using std::vector;

class Wavefunction
{
    public:
        Wavefunction(int NParticles, int NDimensions):
            m_NParticles(NParticles), m_NDimensions(NDimensions)
        {}

        virtual void setParameters( const vector<double> parameters ) = 0;
        virtual double powers( const vector< vector<double> > & positions) = 0;
        virtual double localEnergy( const vector< vector<double> > & positions) = 0;
        //virtual double kinetic() = 0;
        //virtual double potential() = 0;
        //virtual double distance() = 0;
        //virtual void loopParameters() = 0;
        //virtual std::string save() = 0;

    protected:
        int m_NParticles, m_NDimensions;
};

#endif // end of wavefunction base class
