#include "trialwf2.h"
#include <iostream>

using std::cout;
using std::endl;

trialWF2::trialWF2( int NParticles, int NDimensions ):
Wavefunction( NParticles, NDimensions )
{
    if( NParticles != 2 )
    {
        cout << "system defined for 2 particles" << endl;
        exit(1);
    }

}


void trialWF2::setParameters( vector<double> parameters )
{
    if( parameters.size() < 3 )
    {
        cout << "Needs omega, alpha AND beta" << endl;
    }
    m_omega = parameters[0];
    m_alpha = parameters[1];
    m_beta = parameters[2];
}


double trialWF2::sumSquares( const vector< vector<double> > &positions)
{
    double sumSquares = 0;
    for( int particle = 0; particle < m_NParticles; particle ++ )
    {
        for( int dimension = 0; dimension < m_NDimensions; dimension++ )
        {
             sumSquares += positions[particle][dimension]
                          *positions[particle][dimension];
        }
    }
    return sumSquares;
}

double trialWF2::distanceSquared( const vector< vector<double> > &positions )
{
    double difference = 0;
    double distanceSquared = 0;
    for( int dimension = 0; dimension < m_NDimensions; dimension++ )
    {
        difference += positions[0][dimension] - positions[1][dimension];
        distanceSquared += difference*difference;
    }
    return  distanceSquared;
}

double trialWF2::localEnergy( const vector< vector<double> > &positions )
{
    vector<double> energy;
    double r12 = distance(positions);
    double betaFract = 1.0/( (1.0 + m_beta*r12 )*(1.0 + m_beta*r12 ) );
    energy.push_back(
            0.5*m_alpha*m_omega*( 3.0 - m_alpha*m_omega )*sumSquares(positions)
          + 0.5*betaFract*( m_alpha*m_omega*r12
                          - 0.5*betaFract
                          - 2.0/r12
                          + 2*m_beta/(1 + m_beta*r12 ) )
                    );
    energy.push_back(
            0.5*m_omega*m_omega*sumSquares(positions) + 1.0/r12
                    );
    m_kinetic = energy[0];
    m_potential = energy[1];
    m_distance = r12;
    return energy[0]+energy[1];
}
