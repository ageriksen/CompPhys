#include "vmc.h"
#include "storage.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>
#include <mpi.h>

using std::to_string;
using std::cout;
using std::endl;

//------------------------------------------------------
//  VRUN MC
//------------------------------------------------------
void VMC::runVMC( int MCCycles, double steplength )
{
    /*
     * Implementation of the variational Monte Carlo method for finding the
     * optimal trial wavefunction of a given family.
     */

    //------------------------------------------------------
    //initialize ratio variable
    //------------------------------------------------------
    double exponent = 0;
    double ratio = 0;

    //------------------------------------------------------
    //ensure energies to be stored are nulled out before the start.
    //------------------------------------------------------
    double localEnergy = 0;
    m_energy = 0;
    m_energySquared = 0;
    //reset acceptance as well
    m_acceptRatio = 0;
    m_distance = 0;
    m_kinetic = 0;
    m_potential = 0;

    //------------------------------------------------------
    // RNG setup
    //------------------------------------------------------
    std::random_device rd;
    std::mt19937_64 engine(rd() + m_rank);
    std::uniform_real_distribution<double> uniformDistribution(-1, 1);
    std::uniform_real_distribution<double> acceptanceDistribution(0, 1);

    //------------------------------------------------------
    // initial value for system
    //------------------------------------------------------
    for( int particle = 0; particle < m_NParticles; particle ++ )
    {
        for( int dimension = 0; dimension < m_NDimensions; dimension ++ )
        {
            m_positionsOld[particle][dimension]
                = uniformDistribution( engine );
            m_positionsNew[particle][dimension]
                = m_positionsOld[particle][dimension];
        }
    }

    //------------------------------------------------------
    // MCCycles split amongst processors:
    //------------------------------------------------------
    m_MCCycles = MCCycles / m_processors;
    m_MCCyclesFull = MCCycles;

    //------------------------------------------------------
    // initial wavefunction exponent
    //------------------------------------------------------
    m_oldExponent = m_WF -> powers(m_positionsOld);

    //------------------------------------------------------
    //  MC simulation over total cycles
    //------------------------------------------------------
    unsigned long long acceptance = 0;
    for( int cycle = 0; cycle < m_MCCycles; cycle ++ )
    {
        // runs each particle through the metropolis sampling rule
        for( int particle = 0; particle < m_NParticles; particle ++ )
        {
            //Suggest update
            for( int dimension = 0; dimension < m_NDimensions; dimension ++ )
            {
                m_positionsNew[particle][dimension] = m_positionsOld[particle][dimension]
                                                    + steplength*uniformDistribution( engine );
            } // end update suggestion
            // and finds new wavefunction exponent
            m_newExponent = m_WF -> powers(m_positionsNew);

            // updates ratio:
            exponent = 2.0*m_newExponent - 2.0*m_oldExponent;
            ratio = std::exp(exponent);

            if( acceptanceDistribution( engine ) <= ratio )
            {// accepting suggestion
                m_oldExponent = m_newExponent;
                for( int dimension = 0; dimension < m_NDimensions; dimension ++ )
                {
                    m_positionsOld[particle][dimension] = m_positionsNew[particle][dimension];
                }
                // Update acceptance:
                acceptance++;
            }
            else
            {//reset parameters
                m_newExponent = m_oldExponent;
                for( int dimension = 0; dimension < m_NDimensions; dimension ++ )
                {
                    m_positionsNew[particle][dimension] = m_positionsOld[particle][dimension];
                }
            }
        } // end of metropolis per particle

        //------------------------------------------------------
        // New local energy
        //------------------------------------------------------
        localEnergy = m_WF -> localEnergy(m_positionsNew);

        m_energy += localEnergy;
        m_energySquared += localEnergy*localEnergy;
        m_distance += m_WF->distance();
        m_kinetic += m_WF->kinetic();
        m_potential += m_WF->potential();


    } // end of MCCycles
    m_acceptRatio = double(acceptance)/double(m_NParticles*m_MCCycles);


    //------------------------------------------------------
    // Gathering data from processors
    //------------------------------------------------------
    double tmpEnergy = 0;
    double tmpEnergySquared = 0;
    double tmpAcceptRatio = 0;
    double tmpDistance = 0;
    double tmpKinetic = 0;
    double tmpPotential = 0;

    MPI_Reduce( &m_energy,          &tmpEnergy,         1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce( &m_energySquared,   &tmpEnergySquared,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce( &m_acceptRatio,     &tmpAcceptRatio,    1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce( &m_distance,        &tmpDistance,       1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce( &m_kinetic,         &tmpKinetic,        1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce( &m_potential,       &tmpPotential,      1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    m_energy = tmpEnergy;
    m_energySquared = tmpEnergySquared;
    m_acceptRatio = tmpAcceptRatio/m_processors;
    m_distance = tmpDistance/double(m_MCCyclesFull*m_processors);
    m_kinetic = tmpKinetic/double(m_MCCyclesFull*m_processors);
    m_potential= tmpPotential/double(m_MCCyclesFull*m_processors);

    m_variance = ( m_energySquared - m_energy*m_energy/double(m_MCCyclesFull*m_processors))/double(m_MCCyclesFull*m_processors);
    m_energyMean = m_energy/double(m_MCCyclesFull*m_processors);


} // end runVMC


//------------------------------------------------------
//  STEP FINDER
//------------------------------------------------------
double VMC::stepFinder( vector<double> param )
{
    // variable setup
    double trialStep;
    int quickMC = 1e5; // only need acceptance ratio, not exact results
    double acceptanceMin = 0.48; // lower bound on acceptance ratio
    double acceptanceMax = 0.53; // upper bound on acceptance ratio
    double currentAccept;
    double deltaMin = 0.01/double(param[0]); // min delta value tested and step
    double deltaMax = 10/double(param[0]);  // max delta value tested
    std::cout << "min and step: " << deltaMin << ", max: " << deltaMax <<  std::endl;

    //testing steplengths
    for( double delta = deltaMin; delta < deltaMax; delta += deltaMin )
    {
        m_WF -> setParameters( param );
        this -> runVMC( quickMC, delta );

        //std::cout << "delta, acceptratio: " << delta << ", " << m_acceptRatio << "\n";

        if(  m_acceptRatio > acceptanceMin && m_acceptRatio < acceptanceMax)
        {
            currentAccept = m_acceptRatio;
            trialStep = delta;
            std::cout << "steplength: " << trialStep << "\n"
                      << "acceptRatio: " << currentAccept << std::endl;
            clean();
            return trialStep;
        }
    }
    if( m_rank == 0 )
    {
        std::cout << "could not find a good steplength" << std::endl;
//        exit(0);
    }
    clean();
    return 5;
}// end stepFinder

//------------------------------------------------------
// OPTIMIZE
//------------------------------------------------------
vector<double> VMC::optimize
(
    vector< vector<double> > parameters,
    vector<int> range,
    double delta,
    double alpha0
)
{
    m_stepLength = delta;
    m_parameters = parameters;
    vector<double> minima(2);
    minima[0] = alpha0;
    minima[1] = beta0( minima[0], m_MCCyclesFull );
    for( unsigned int i = 1; i < range.size(); i++ )
    {
        double alpha = findAlpha( minima );
        double beta = findBeta( minima );
        minima[0] = alpha;
        minima[1] = beta;
    }
    return minima;
} // end optimize

//------------------------------------------------------
// FINDERS
//------------------------------------------------------
double VMC::alpha0
(
    vector< vector<double> > parameters,
    int MCCycles,
    string baseName
)
{ // no change in minimal alpha per omega. Set it to 1.0

    cout << "finding alpha0" << endl;

    vector<double> param(2);
    param[0] = parameters[0].back(); // omega == 1.0
    string tmpLine = "";

    double delta = parameters[1].back();
    double energyMin = m_energyMean;
    double alpha0 = 0;

    for( double alpha: parameters[2] )
    {
        param[1] = alpha;

        m_WF -> setParameters( param );
        runVMC( MCCycles, delta );

        if( m_energyMean < energyMin )
        {
            alpha0 = param[1];
        }

    }

if( m_rank == 0 )
    {
        tmpLine = to_string(alpha0) + " "
                + to_string(m_energyMean) + " "
                + to_string(m_variance) + " "
                + to_string(m_acceptRatio);

        storage variableSaver(baseName + "alpha0.dat");
        variableSaver.out();
        variableSaver.dat(tmpLine);
        cout << tmpLine << endl;
        variableSaver.close();

        cout << "-------------------------------------" << endl;
        cout << "saving alpha0" << alpha0 << endl;
        storage alpha0Saver("./resources/alpha0.dat");
        alpha0Saver.out();
        alpha0Saver.dat(to_string(alpha0));
        alpha0Saver.close();
    }
    return alpha0;
} // end alpha0


double VMC::beta0
(
    double alpha0,
    int MCCycles
)
{ // no change in minimal alpha per omega. Set it to 1.0

    cout << "finding beta0" << endl;

    vector<double> param(3);
    param[0] = m_parameters[0].back(); // omega == 1.0
    param[1] = alpha0;

    double delta = m_parameters[1].back();
    double energyMin = m_energyMean;
    double beta0 = 0;

    string tmpLine = "";

    for( double beta: m_parameters[3] )
    {
        param[3] = beta;
        m_WF -> setParameters( param );
        runVMC( MCCycles, delta );

        if( m_variance < energyMin )
        {
            beta0 = param[3];
        }

    }

    cout << "beta0 = " << beta0 << endl;

    return beta0;
} // end alpha0

double VMC::findAlpha( vector<double> param)
{
    //double varianceMin;
    //double alphaMin;
    //vector<double> param(2);
    //param[1] = beta;
    //for( double alpha: m_parameters[1] )
    //{
    //    param[0] = alpha;
    //    m_WF -> setParameters( param );
    //    runVMC( m_MCCyclesFull, m_stepLength );
    //    if( m_variance < varianceMin )
    //    {
    //        varianceMin = m_variance;
    //        alphaMin = param[0];
    //    }
    //}

    double oldMin = param[0];
    int counter=0;
    double epsilon = 1;
    double ratio=epsilon*2;
    double energyOld = 0;
    double width = param[0]*0.05; // half the width of integration;
    double delta = width*0.01;
    param[0] -= width;

    cout << "entering alpha whileLoop" << endl;
    while( !(ratio<epsilon) )
    {
        if( counter > 100 )
        {
            cout << "Over count limit. terminating alpha: "
                 << param[0] << ". resetting to " << oldMin << endl;
            param[0] = oldMin;
            break;
        }
        energyOld = m_energy;
        oldMin = param[0];
        param[0] += delta;
        m_WF -> setParameters( param );
        runVMC( m_MCCyclesFull, m_stepLength );
        ratio = m_energy/energyOld;
        if(ratio > 1)
        {
            cout << "turning back" << endl;
            delta *= -0.5;
            param[0] += delta;
            m_WF -> setParameters( param );
            runVMC( m_MCCyclesFull, m_stepLength );
            ratio = m_energy/energyOld;
        }
        counter++;
    }
    cout << "shift in alpha: " << oldMin - param[0] << endl;
    return param[0];
} // end findAlpha


double VMC::findBeta( vector<double> param )
{
    //double varianceMin;
    //double betaMin;
    //vector<double> param(2);
    //param[1] = alpha;
    //for( double beta: m_parameters[2] )
    //{
    //    param[1] = beta;
    //    m_WF -> setParameters( param );
    //    runVMC( m_MCCyclesFull, m_stepLength );
    //    if( m_variance < varianceMin )
    //    {
    //        varianceMin = m_variance;
    //        betaMin = param[1];
    //    }
    //}
    double firstMin = param[1];
    double oldMin = param[1];
    int counter=0;
    double epsilon = 1;
    double ratio=epsilon*2;
    double energyOld = 0;
    double width = param[1]*0.05; // half the width of integration;
    double delta = width*0.01;
    param[1] -= width;

    cout << "entering beta whileLoop" << endl;
    while( !(ratio<epsilon) )
    {
        if( counter > 100 )
        {
            cout << "Over count limit. terminating at beta: "
                 << param[1] << ". resetting to " << oldMin << endl;
            param[1] = oldMin;
            break;
        }
        energyOld = m_energy;
        oldMin = param[1];
        param[1] += delta;
        m_WF -> setParameters( param );
        runVMC( m_MCCyclesFull, m_stepLength );
        ratio = m_energy/energyOld;
        if(ratio > 1)
        {
            cout << "turning back" << endl;
            delta *= -0.5;
            param[1] += delta;
            m_WF -> setParameters( param );
            runVMC( m_MCCyclesFull, m_stepLength );
            ratio = m_energy/energyOld;
        }
        counter++;
    }
    cout << "shift in beta: " << firstMin - param[1] << endl;
    return param[1];
} // end findBeta
