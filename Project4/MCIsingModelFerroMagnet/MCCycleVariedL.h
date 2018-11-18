#ifndef MCCYCLE_H
#define MCCYCLE_H

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
        string,

        int ,
        int ,
        int ,

        double ,

        vec & ,
        vec & ,
        vec & ,
        vec & ,
        vec & ,
        vec & ,
        vec & ,
        vec & ,
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
        mat ,
        int
        );

#endif
