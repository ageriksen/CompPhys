#include "trialwf1full.h"
#include <iostream>

trialWF1Full::trialWF1Full
(
    int NParticles,
    int NDimensions
) : Wavefunction( NParticles, NDimensions )
{
    if( NParticles != 2 )
    {
        std::cout
            << "system defined for 2 particles"
            << std::endl;
        exit(1);
    }
} // end constructor

double trialWF1Full::sumSquares
(
 const vector< vector<double> > & positions
)
{
    double sumSquares = 0;

    for( int particle = 0; particle < m_NParticles; particle++ )
    {
        for( int dimension = 0; dimension < m_NDimensions; dimension++)
        {
            sumSquares += positions[particle][dimension]
                         *positions[particle][dimension];
        }
    }
    return sumSquares;
}

double trialWF1Full::distanceSquared(const vector<vector<double> > &positions)
{
    double difference;
    for( int dimension = 0; dimension < m_NDimensions; dimension++ )
    {//possibly implement 2 forloops, 1 with i from 0-NParticles,
    // and 1 j from i+1-NParticles. This is redundant in case w/
    // restrictions to merely 2 particles.
        difference = positions[0][dimension] - positions[1][dimension];
        m_distanceSquared += difference*difference;
    }
    return m_distanceSquared;
}
