#include "MCCycle.h"
#include <iostream>

inline int PeriodicBoundary(int i, int limit, int add)
{
    return (i + limit + add) % (limit);
}

int main(int argc, char *argv[]){

    string runName;
    int NSpins, MCCycles;
    double initialTemp, finalTemp, TempStep;

    cout << "please input: \n"
            << "read runName, Number-of-Spins, MC-cycles, "
            << "initial-temperature, final-temperature "
            << "and temperature-step." << endl;
    cin >> runName;
    cin >> NSpins;
    cin >> MCCycles;
    cin >> initialTemp;
    cin >> finalTemp;
    cin >> TempStep;

    // run MC cycles over temperature range
    for( double Temp = initialTemp; Temp <= finalTemp; Temp += TempStep )
    {
        string Filename = runName+std::to_string(Temp)+".bin";
        vec Expectationvalues = zeros<mat>(6);
        MonteCarloMetropolis( runName, NSpins, MCCycles, Temp, Expectationvalues );
//        Expectationvalues /= MCCycles;
//        fileDump( Filename, Expectationvalues, Expectationvalues.size() );
    }
    return 0;
} // end main

void MonteCarloMetropolis(
        string runName,
        int NSpins,
        int MCCycles,
        double Temp,
        vec & Expectationvalues)
{
    // Seed and Mersenne
    std::random_device rd;
    std::mt19937_64 gen(rd());
    // Uniform distribution setup 1 >= x >= 0
    std::uniform_real_distribution<double> RNG(0.0, 1.0);

    //Lattice initialization
    mat SpinMatrix = zeros<mat>(NSpins, NSpins);
    //Energy and magnetization
    double Energy = 0.; double  MagneticMoment = 0.;
    //Expectationvalues array
    initializeLattice( NSpins, SpinMatrix, Energy, MagneticMoment);
    //Array for possible changes in energy
    vec EnergyProb = zeros<mat>(17);
    for( int dE = -8; dE <= 8; dE += 4 ) EnergyProb(dE+8) = exp(-dE/Temp);
    int tenth = MCCycles*0.1;

    //Begin Monte Carlo cycle
    for( int cycle = 1; cycle <= MCCycles; cycle ++)
    {
        //sweep lattice
        for( int count = 0; count <= (NSpins*NSpins); count ++ )
        {
            int x = (int)( RNG(gen)*(double)NSpins );
            int y = (int)( RNG(gen)*(double)NSpins );
            int DeltaE = 2*SpinMatrix(x, y)*
                (
                 SpinMatrix( x, PeriodicBoundary( y, NSpins, -1) )
                 + SpinMatrix( PeriodicBoundary( x, NSpins, -1), y )
                 + SpinMatrix( x, PeriodicBoundary( y, NSpins, 1) )
                 + SpinMatrix( PeriodicBoundary( y, NSpins, 1), y )
                );
            if( RNG(gen) <= EnergyProb(DeltaE + 8) )
            {
                SpinMatrix(x, y) *= -1.; // flip and accept spin
                MagneticMoment += (double) 2*SpinMatrix(x, y);
                Energy += (double)DeltaE;
            }
        }
        //Update expectation values, local node
        Expectationvalues(0) = cycle;
        Expectationvalues(1) += Energy;
        Expectationvalues(2) += Energy*Energy;
        Expectationvalues(3) += MagneticMoment;
        Expectationvalues(4) += MagneticMoment*MagneticMoment;
        Expectationvalues(5) += fabs(MagneticMoment);

        if( cycle % tenth == 0)
        {
            fileDump( runName+"sweeps"+std::to_string(cycle)+".bin", Expectationvalues/(cycle*NSpins*NSpins), Expectationvalues.size() );
        }
    }

} // end MonteCarloMetropolis

void initializeLattice( int NSpins, mat & SpinMatrix, double & Energy, double & MagneticMoment )
{
    /*
     *Implement initialization and starting configuration with call to RNG function in future.
     */

    // initialize spin matrix and magnetization, ordered
    for( int x = 0; x < NSpins; x ++ )
    {
        for( int y = 0; y < NSpins; y ++ )
        {
            SpinMatrix(x, y) = 1.; //g.s. spin orientation
            MagneticMoment += (double)SpinMatrix(x, y);
        }
    }
    // initial energy
    for( int x = 0; x < NSpins; x ++ )
    {
        for( int y = 0; y < NSpins; y ++ )
        {
            Energy -= (double)SpinMatrix(x,y)*
                (
                 SpinMatrix( PeriodicBoundary( x, NSpins, -1 ), y )
                 + SpinMatrix( x, PeriodicBoundary( y, NSpins, -1 ) )
                );
        }
    }
} // end initializeLattice

void fileDump(
        string Filename,
        mat Array,
        int ArrayDim
        )
{
    /*
     * I want an initializer to take a _Filename
     * I want a call open() which opens the file and a close() for closing it
     * I want a dump() that dumps a bitmap
     *
     * First off, dump and make rest as a class
     */

    std::ofstream file( Filename, std::ofstream::binary);

    double *tmp_arr = new double[ArrayDim];
    for( int i = 0; i < ArrayDim; i++)
    {
        tmp_arr[i] = Array[i];
    }

    file.write
        (
            reinterpret_cast<const char*> (tmp_arr), ArrayDim*sizeof(double)
        );

    file.close();
    delete [] tmp_arr;

//    FILE* outfile = fopen(Filename.c_str(), "wb");
//    fwrite( & Array[0], sizeof(double), ArrayDim*ArrayDim, outfile );
//    fclose(outfile);

} // end fileDump
