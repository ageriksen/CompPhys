#ifndef TRIALWF2_H
#define TRIALWF2_H

#include "wavefunction.h"

class trialWF2: public Wavefunction
{
    public:
        trialWF2( int NParticles, int NDimensions);

        void setParameters( vector<double> parameters )
        {
            m_omega = parameters[0];
            m_alpha = parameters[1];
            m_beta = parameters[2];
        }

        double sumSquares( const vector< vector<double> > &positions );
        double distanceSquared( const vector< vector<double> > &positions );
        double distance( const vector< vector<double> > &positions );

        double powers( const vector< vector<double> > &positions );
        double localEnergy( const vector< vector<double> > &positions );

        //getters, setters
        double omega() { return m_omega; }
        double alpha() { return m_alpha; }
        double beta() { return m_beta; }
        double distance() { return m_distance; }
        double kinetic() { return m_kinetic; }
        double potential() { return m_potential; }

    private:
        double m_omega = 0;
        double m_alpha = 0;
        double m_beta = 0;
        double m_distance = 0;
        double m_kinetic = 0;
        double m_potential = 0;
}; // end of trialWF2

inline double trialWF2::distance( const vector< vector<double> > &positions )
{ return std::sqrt( distanceSquared(positions) ); }

inline double trialWF2::powers( const vector< vector<double> > &positions )
{
    return 0.5*( 1.0/( 1.0/distance(positions) + m_beta )
         - m_alpha*m_omega*sumSquares(positions) );
}

#endif // end TRIALWF2_H
