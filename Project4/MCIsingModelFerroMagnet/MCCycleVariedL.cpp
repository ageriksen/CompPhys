#include "MCCycle.h"
#include <iostream>

inline int PeriodicBoundary(int i, int limit, int add)
{
    return (i + limit + add) % (limit);
}

int main(int argc, char *argv[]){

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    string path, runName, mode;
    int NSpins, MCCycles, Equilibrium;
    double initialTemp, finalTemp, TempStep;//, Rate;
    vec ExpectEnergy, ExpectEnergySquared, ExpectMagnet, ExpectMagnetSquared, AbsMagnet;
    vec Accepted, MCTime, EnergyArray;
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cin >> path;
    cin >> runName;
    cin >> mode;
    cin >> NSpins;
    cin >> MCCycles;
    cin >> Equilibrium;
    //cin >> Rate;
    cin >> initialTemp;
    cin >> finalTemp;
    cin >> TempStep;
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // { Rate of stotrage indefinitely suspended}
    EnergyArray = zeros<vec>(MCCycles - Equilibrium);
    Accepted = zeros<vec>(MCCycles - Equilibrium);
    ExpectEnergy = zeros<vec>(MCCycles - Equilibrium);
    ExpectEnergySquared = zeros<vec>(MCCycles - Equilibrium);
    ExpectMagnet = zeros<vec>(MCCycles - Equilibrium);
    ExpectMagnetSquared = zeros<vec>(MCCycles - Equilibrium);
    AbsMagnet = zeros<vec>(MCCycles - Equilibrium);
    MCTime = zeros<vec>(MCCycles - Equilibrium);
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // run MC cycles over temperature range
    int TempCount = 0;
    for( double Temp = initialTemp; Temp <= finalTemp; Temp += TempStep )
    {
        cout << "T : " << Temp << "\n";
        TempCount++;
        vec Expectationvalues = zeros<mat>(5);
        MonteCarloMetropolis(
                runName,
                NSpins,
                MCCycles,
                Equilibrium,
                Temp,
                mode,
                EnergyArray,
                Expectationvalues,
                ExpectEnergy,
                ExpectEnergySquared,
                ExpectMagnet,
                ExpectMagnetSquared,
                AbsMagnet,
                Accepted,
                MCTime
               );

        // printing results to standard out and writing arrays to file:
        Expectationvalues /= (MCCycles);

        cout << " <E>/(N^2) : " << Expectationvalues(0)/(NSpins*NSpins) << "\n"
             << " <|M|>/(N^2) : " << Expectationvalues(4)/(NSpins*NSpins) << "\n"
             << " Cv*(k*T*T)/(N^2) : " << (Expectationvalues(1) - Expectationvalues(0)*Expectationvalues(0))/(Temp*Temp*NSpins*NSpins) << "\n"
             << " X*(k*T)/(N^2) : " << (Expectationvalues(3) - Expectationvalues(4)*Expectationvalues(4))/(Temp*NSpins*NSpins) << endl;

        cout << "Final expectation values, separate: \n"
             << Expectationvalues
        //setting up strings to reduce errors in writing:
        string TempString = std::to_string(TempCount);
        string Filename = "ExpectationValuesTemp"+TempString+runName;
        string FinalFilename = path+"ExpectationValuesTemp"+TempString+"Final"+runName;
        // writing to file:
        fileDump( path+"Energy"+Filename, ExpectEnergy, ExpectEnergy.size() );
        fileDump( path+"EnergySquared"+Filename, ExpectEnergySquared, ExpectEnergySquared.size() );
        fileDump( path+"Magnetization"+Filename, ExpectMagnet, ExpectMagnet.size() );
        fileDump( path+"MagnetizationSquared"+Filename, ExpectMagnetSquared, ExpectMagnetSquared.size() );
        fileDump( path+"MagnetizationAbsoluteValue"+Filename, AbsMagnet, AbsMagnet.size() );
        fileDump( path+"AcceptedConfigurationsTemp"+TempString+runName, Accepted, Accepted.size());

        fileDump( FinalFilename, Expectationvalues, Expectationvalues.size() );
    }
    return 0;
} // end main
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void MonteCarloMetropolis(
        string runName,
        int NSpins,
        int MCCycles,
        int Equilibrium,
        double Temp,
        //double Rate,
        string mode,
        vec & EnergyArray,
        vec & Expectationvalues,
        vec & ExpectEnergy,
        vec & ExpectEnergySquared,
        vec & ExpectMagnet,
        vec & ExpectMagnetSquared,
        vec & AbsMagnet,
        vec & Accepted,
        vec & MCTime
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
    //
    //
    ////Rate of storage {INDEFINITE SUSPENTION}
    //int step = (int)(Rate*(double)MCCycles);

    //Begin Monte Carlo cycle
    for( int cycle = 0; cycle < MCCycles; cycle ++)
    {
        int counter = 0;
        //sweep lattice
        for( int count = 0; count < (NSpins*NSpins); count ++ )
        {
            int x = (int)( RNG(gen)*(double)NSpins );
            int y = (int)( RNG(gen)*(double)NSpins );

            int DeltaE = 2*Lattice(x, y)*
                (
                 Lattice( x, PeriodicBoundary( y, NSpins, -1) )
                 + Lattice( PeriodicBoundary( x, NSpins, -1), y )
                 + Lattice( x, PeriodicBoundary( y, NSpins, 1) )
                 + Lattice( PeriodicBoundary( x,  NSpins, 1), y )
                );
            if( RNG(gen) <= EnergyProb(DeltaE + 8) )
            {
                Lattice(x, y) *= -1.; // flip and accept spin
                MagneticMoment += (double) 2*Lattice(x, y);
                Energy += (double)DeltaE;
                counter += 1;
            }
        }
        //Update expectation values, local node
        Expectationvalues(0) += Energy;
        Expectationvalues(1) += Energy*Energy;
        Expectationvalues(2) += MagneticMoment;
        Expectationvalues(3) += MagneticMoment*MagneticMoment;
        Expectationvalues(4) += fabs(MagneticMoment);

        // For post exercise C, find ~equilibration time, then start measuring after.
        if( cycle >= Equilibrium )
        {
            int index = cycle - Equilibrium;
            // {RATE OF STORAGE INDEFINITELY SUSPENDED}

            EnergyArray(index) = Energy;
            Accepted(index) = counter;
            double denominator = (cycle+1)*NSpins*NSpins;

            ExpectEnergy(index) = Expectationvalues(0)/denominator;
            ExpectEnergySquared(index) = Expectationvalues(1)/denominator;
            ExpectMagnet(index) = Expectationvalues(2)/denominator;
            ExpectMagnetSquared(index) = Expectationvalues(3)/denominator;
            AbsMagnet(index) = Expectationvalues(4)/denominator;
        }
    }

} // end MonteCarloMetropolis
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void initializeLattice( int NSpins, mat & Lattice, double & Energy, double & MagneticMoment, string mode)
{
    /*
     * Implement initialization and starting configuration with call to RNG function in future.
     */

    if (strcmp(mode.c_str(),"Up") == 0)
    {   // Fill lattice with spin up
        Lattice.ones();
        MagneticMoment = Lattice.size();
    }
    else if (strcmp(mode.c_str(),"Down") == 0)
    {   // Fill lattice with spin down
        Lattice.ones();
        Lattice = Lattice*-1;
        MagneticMoment = -1.0*Lattice.size();
    }
    else if (strcmp(mode.c_str(),"Random") == 0)
    {   // Fill lattice with random spin directions
        mat randomLattice(size(Lattice),fill::randu); randomLattice = sign(randomLattice*2-1);
        Lattice = randomLattice;
        MagneticMoment = 0.0;
        for (int x = 0; x < NSpins; ++x)
        {
            for (int y = 0; y < NSpins; ++y)
            {
                MagneticMoment += Lattice(x,y);
            }
        }
    }
    else
    {
        cout << "Did you remember proper capitalization? Either 'Up', 'Down' or 'Random'" << endl;
    }

    Energy = 0;
    for (int i = 0; i< NSpins; i++)
    {
        for (int j = 0; j< NSpins; j++)
        {
            Energy -= (double)(Lattice(i,j)*(
                        Lattice(PeriodicBoundary(i, NSpins, 1), j)
                      + Lattice(i, PeriodicBoundary(j, NSpins, 1))
                      ));
        }
    }
    cout << "initial energy per particle is: " << Energy/(NSpins*NSpins) << endl;
} // end initializeLattice
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

} // end fileDump
