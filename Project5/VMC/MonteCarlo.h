#ifndef MONTECARLO_H
#define MONTECARLO_H

#include "Random.h"

class MonteCarlo
{
    public:
        MonteCarlo();  // constructor
        ~MonteCarlo(); // destructor

        void sweep(); // 1 sweep over system

        void cycle(); // cycle N sweeps

    private:
}; // end MonteCarlo class

#endif // end MonteCarlo header
