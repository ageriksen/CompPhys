#include "MCCycle.h"
#include <iostream>

inline int PeriodicBoundary(int i, int limit, int add)
{
    return (i + limit + add) % (limit);
}

int main(int argc, char *argv[]){

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    string path, runName, mode;
    int NSpins, MCCycles;
    double initialTemp, finalTemp, TempStep, Rate;
    vec ExpectEnergy, ExpectEnergySquared, ExpectMagnet, ExpectMagnetSquared, AbsMagnet;
    vec Accepted, MCTime;
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cout << "please input: \n"
            << "path, runName, mode, Nspins, MCCycles, Rate, InitialTemp, "
            << "finalTemp and TempStep." << endl;
    cin >> path;
    cin >> runName;
    cin >> mode;
    cin >> NSpins;
    cin >> MCCycles;
    cin >> Rate;
    cin >> initialTemp;
    cin >> finalTemp;
    cin >> TempStep;
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cout << "necessary values read." << endl;

    cout << "with a fraction of " << Rate << endl;
    cout << "and a max number of cycles " << MCCycles << endl;
    cout << "\n \n The length of our arrays are: " << MCCycles*Rate << "\n" << endl;
    // Initializing arrays for writing:
    Accepted = zeros<vec>(MCCycles*Rate);
    ExpectEnergy = zeros<vec>(MCCycles*Rate);
    ExpectEnergySquared = zeros<vec>(MCCycles*Rate);
    ExpectMagnet = zeros<vec>(MCCycles*Rate);
    ExpectMagnetSquared = zeros<vec>(MCCycles*Rate);
    AbsMagnet = zeros<vec>(MCCycles*Rate);
    MCTime = zeros<vec>(MCCycles*Rate);
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cout << "arrays initialized. " << endl;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // run MC cycles over temperature range
    int TempCount = 0;
    for( double Temp = initialTemp; Temp <= finalTemp; Temp += TempStep )
    {
        cout << "let's run the Monte Carlo cycles." << endl;
        TempCount++;
        vec Expectationvalues = zeros<mat>(5);
        MonteCarloMetropolis(
                runName,
                NSpins,
                MCCycles,
                Rate,
                Temp,
                mode,
                Expectationvalues,
                ExpectEnergy,
                ExpectEnergySquared,
                ExpectMagnet,
                ExpectMagnetSquared,
                AbsMagnet,
                Accepted,
                MCTime
               );

        cout << "Energy expectation value Array: \n"
             << ExpectEnergy << endl;
        // printing results to standard out and writing arrays to file:
        Expectationvalues /= (MCCycles);
        cout << "Final expexted energy per spin: " << Expectationvalues(0)/(NSpins*NSpins)
             << " at temperature " << Temp*initialTemp << endl;
        cout << "Final expectation values, separate: \n"
             << Expectationvalues
             << "Energy per spin: " << Expectationvalues(0)/(NSpins*NSpins) << "\n"
             << "mean abosolute magnetization per spin: " << Expectationvalues(4)/(NSpins*NSpins) << "\n"
             << "heat capacity Cv*(kT*T) per spin: " << (Expectationvalues(1) - Expectationvalues(0)*Expectationvalues(0))/(NSpins*NSpins) << "\n"
             << "Susceptibility Chi*kT per spin: " << (Expectationvalues(3) - Expectationvalues(2)*Expectationvalues(2))/(NSpins*NSpins) << endl;
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
        double Temp,
        double Rate,
        string mode,
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
    //Rate of storage
//    int Rate = (int)((float)MCCycles*(float)Rate);
    int step = (int)(Rate*(double)MCCycles);
    cout << "we want to save every " << step << "step" << endl;

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

//        cout << "we're on cycle: " << cycle << "\n"
//             << "and we want to check the rmainder of cycle % step with a rate: " << step << ": \n"
//             << cycle % step << "\n";
        if( cycle % step == 0 )
        {
            cout << "and so we save on cycle "<< cycle << endl;
            Accepted(cycle) = counter;
            double denominator = (cycle+1)*NSpins*NSpins;

            ExpectEnergy(cycle) = Expectationvalues(0)/denominator;
            ExpectEnergySquared(cycle) = Expectationvalues(1)/denominator;
            ExpectMagnet(cycle) = Expectationvalues(2)/denominator;
            ExpectMagnetSquared(cycle) = Expectationvalues(3)/denominator;
            AbsMagnet(cycle) = Expectationvalues(4)/denominator;
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
