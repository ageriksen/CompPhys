#include "MCCycleEnergies.h"
#include <iostream>

inline int PeriodicBoundary(int i, int limit, int add)
{
    return (i + limit + add) % (limit);
}

int main(int argc, char *argv[]){

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    string path, mode;
    //
    int Lmin, Lmax, Lstep;
    int MCCycles, Equilibrium;
    //
    double initialTemp, finalTemp, TempStep;
    //
    vec Energies;
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // strings:
    cin >> path;
    cin >> mode;
    // integers:
    cin >> Lmin;
    cin >> Lmax;
    cin >> Lstep;
    cin >> MCCycles;
    cin >> Equilibrium;
    // floats:
    cin >> initialTemp;
    cin >> finalTemp;
    cin >> TempStep;
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clock_t timeStart, timeFinish;
    double timeused;
    timeStart = clock();
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for( int NSpins = Lmin; NSpins <= Lmax; NSpins += Lstep )
    {
        //
        Energies = zeros<vec>(MCCycles - Equilibrium);
        //
        cout << "==================================\n"
             << "L : " << NSpins << "\n";
        //

        // run MC cycles over temperature range
        int TempCount = 0;
        for( double Temp = initialTemp; Temp <= finalTemp; Temp += TempStep )
        {
            //
            cout << "T : " << Temp << "\n";
            //
            vec Expectationvalues = zeros<mat>(4);
            //
            MonteCarloMetropolis(
                    mode,
                    //
                    NSpins,
                    MCCycles,
                    Equilibrium,
                    //
                    Temp,
                    //
                    Expectationvalues,
                    Energies
                   );
            // printing results to standard out and writing arrays to file:
            Expectationvalues /= (MCCycles);
            double VarianceE = (
                      Expectationvalues(1)
                    - (  Expectationvalues(0)*Expectationvalues(0) )
                    );
            //
            cout << "| Variance, Energy : " << VarianceE << "\n";
            //
            //setting up strings to reduce errors in writing:
            //
            //Making temp and L more print friendly:
            //
            string Lstring;
            if( NSpins < 10 )
            {
                Lstring = "0"+std::to_string(NSpins);
            }
            else
            {
                Lstring = std::to_string(NSpins);
            }
            //
            //string TempString = std::to_string(int(round(
            //        TempCount*TempStep + initialTemp
            //        )));
            string TempString = std::to_string( int(round( Temp*10 )) );
            string Fname = mode+"L"+Lstring+"T"+TempString+".bin";
            // writing to file:
            fileDump( path+"Energies"+Fname, Energies );
            //
            //
            timeFinish = clock();
            timeused = (double)( timeFinish - timeStart )/CLOCKS_PER_SEC;
            cout << " time: " << timeused << endl;
        }
            //
            TempCount++;
    }// End, spins states
    return 0;
} // end main

/*
 * Following are 3 functions:
 * MonteCarloMetropolis,
 * initializeLattice and
 * fileDump.
 */
//
//
//

void MonteCarloMetropolis(
        string mode,
        //
        int NSpins,
        int MCCycles,
        int Equilibrium,
        //
        double Temp,
        //
        vec & Expectationvalues,
        vec & Energies
        )
{
    //
    // Seed and Mersenne
    std::random_device rd;
    std::mt19937_64 gen(rd());
    //
    // Uniform distribution setup 1 >= x >= 0
    std::uniform_real_distribution<double> RNG(0.0, 1.0);

    //Lattice initialization
    mat Lattice = zeros<mat>(NSpins, NSpins);
    //Energy and magnetization
    double Energy = 0.; double  MagneticMoment = 0.;
    //
    //Expectationvalues array
    initializeLattice( NSpins, Lattice, Energy, MagneticMoment, mode);
    //
    //Array for possible changes in energy
    vec EnergyProb = zeros<mat>(17);
    for( int dE = -8; dE <= 8; dE += 4 ) EnergyProb(dE+8) = exp(-dE/Temp);

    //Begin Monte Carlo cycle'
    int x, y, DeltaE;
    int Accepted = 0;
    for( int cycle = 0; cycle < Equilibrium; cycle ++)
    { // Equilibration loop, to loop over the necessary:
        //sweep lattice
        for( int count = 0; count < (NSpins*NSpins); count ++ )
        {
            x = (int)( RNG(gen)*(double)NSpins );
            y = (int)( RNG(gen)*(double)NSpins );

            DeltaE = 2*Lattice(x, y)*
                (
                   Lattice( x, PeriodicBoundary( y, NSpins, -1) )
                 + Lattice( PeriodicBoundary( x, NSpins, -1), y )
                 + Lattice( x, PeriodicBoundary( y, NSpins,  1) )
                 + Lattice( PeriodicBoundary( x,  NSpins, 1), y )
                );
            if( RNG(gen) <= EnergyProb(DeltaE + 8) )
            {
                Lattice(x, y) *= -1.; // flip and accept spin
                Energy += (double)DeltaE;
            }
        }
    } // end equilibration loop

    for( int cycle = Equilibrium; cycle < MCCycles; cycle ++)
    {
        //sweep lattice
        for( int count = 0; count < (NSpins*NSpins); count ++ )
        {
            x = (int)( RNG(gen)*(double)NSpins );
            y = (int)( RNG(gen)*(double)NSpins );

            DeltaE = 2*Lattice(x, y)*
                (
                   Lattice( x, PeriodicBoundary( y, NSpins, -1) )
                 + Lattice( PeriodicBoundary( x, NSpins, -1), y )
                 + Lattice( x, PeriodicBoundary( y, NSpins,  1) )
                 + Lattice( PeriodicBoundary( x,  NSpins, 1), y )
                );
            if( RNG(gen) <= EnergyProb(DeltaE + 8) )
            {
                Lattice(x, y) *= -1.; // flip and accept spin
                MagneticMoment += (double) 2*Lattice(x, y);
                Energy += (double)DeltaE;
                Accepted ++;
            }
        }
        //Update expectation values, local node
        Expectationvalues.at(0) += Energy;
        Expectationvalues.at(1) += Energy*Energy;
        Expectationvalues.at(2) += MagneticMoment*MagneticMoment;
        Expectationvalues.at(3) += fabs(MagneticMoment);

        int EquilibratedCycle = cycle-Equilibrium;
        Energies(EquilibratedCycle) = Energy;

    }

} // end MonteCarloMetropolis
//
//
//

void initializeLattice( int NSpins, mat & Lattice, double & Energy, double & MagneticMoment, string mode)
{
    /*
     * Initialize Lattice was adapted from Jacod Lilleborg's project code with his permission.
     * Uncertain exactly which folder, but should be in there somewhere:
     * link: https://github.com/Lilleborg/FYS3150-Computational-physics/tree/master/Projects/Project4
     */
//
    if (strcmp(mode.c_str(),"Up") == 0)
    {   // Fill lattice with spin up
        Lattice.ones();
        //
        MagneticMoment = Lattice.size();
    }
    //
    //
    else if (strcmp(mode.c_str(),"Down") == 0)
    {   // Fill lattice with spin down
        Lattice.ones();
        Lattice *= -1;
        //
        MagneticMoment = -1.0*Lattice.size();
    }
    //
    //
    else if (strcmp(mode.c_str(),"Random") == 0)
    {   // Fill lattice with random spin directions
        mat randomLattice(size(Lattice),fill::randu); randomLattice = sign(randomLattice*2-1);
        Lattice = randomLattice;
        //
        MagneticMoment = 0.0;
        for (int x = 0; x < NSpins; ++x)
        {
            for (int y = 0; y < NSpins; ++y)
            {
                MagneticMoment += Lattice(x,y);
            }
        }
    }
    //
    //
    else
    {
        cout
            << "Did you remember proper capitalization? Either 'Up', 'Down' or 'Random'"
            << endl;
    }

    Energy = 0;
    for (int i = 0; i< NSpins; i++)
    {
        for (int j = 0; j< NSpins; j++)
        {
            //
            //
            Energy -= (double)(Lattice(i,j)*(
                        Lattice(PeriodicBoundary(i, NSpins, 1), j)
                      + Lattice(i, PeriodicBoundary(j, NSpins, 1))
                      ));
            //
            //
        }
    }
//    cout << "E0/(N^2) : " << Energy/(NSpins*NSpins) << endl;
} // end initializeLattice
//
//
//

void fileDump(
        string Filename,
        mat Array
        )
{
    /*
     * I want an initializer to take a _Filename
     * I want a call open() which opens the file and a close() for closing it
     * I want a dump() that dumps a bitmap
     *
     * First off, dump and make rest as a class
     */
//
    std::ofstream file( Filename, std::ofstream::binary);
    //
    int ArrayDim = Array.size();
    double *tmp_arr = new double[ArrayDim];
    for( int i = 0; i < ArrayDim; i++)
    {
        tmp_arr[i] = Array[i];
    }
    //
    //
    file.write
        (
            reinterpret_cast<const char*> (tmp_arr), ArrayDim*sizeof(double)
        );
    //
    //
    file.close();
    delete [] tmp_arr;
//
} // end fileDump
