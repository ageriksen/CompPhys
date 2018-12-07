#include "twoparticlenoninteractingwf.h"
#include <iostream>

TwoParticleNonInteractingWF::TwoParticleNonInteractingWF
(
    int NParticles,
    int NDimensions
):
    Wavefunction( NParticles, NDimensions )
{
    if( NParticles != 2 )
    {
        std::cout
            << "system defined for 2 particles"
            << std::endl;
        exit(1);
    }
}

double TwoParticleNonInteractingWF::sumSquares
(
 const arma::Mat<double> & positions
)
{
    double sumSquares = 0;

    for( int particle = 0; particle < m_NParticles; particle++ )
    {
        for( int dimension = 0; dimension < m_NDimensions; dimension++)
        {
            sumSquares += positions(particle, dimension)
                         *positions(particle, dimension);
        }
    }
    return sumSquares;
}
