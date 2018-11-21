#ifndef MCCYCLEENERGIES_H
#define MCCYCLEENERGIES_H

#include <cmath>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <random>
#include <armadillo>
#include <string>

using namespace std;
using namespace arma;

int PeriodicBoundary(
        int ,
        int ,
        int
        );

void MonteCarloMetropolis(
        string ,

        int ,
        int ,
        int ,

        double ,

        vec &,
        vec &
        );

void initializeLattice(
        int ,
        mat & ,
        double & ,
        double & ,
        string
        );

void fileDump(
        string ,
        mat
        );

#endif
