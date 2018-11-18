#include <iostream>
#include <cmath>
#include <cstdlib>
#include "solarsystem.h"
#include "solver.h"
using namespace std;

int main(int numArguments, char **arguments)
{
    int numTimesteps = 1e8;
    int timelength = 100;
    int time_normalization = 365; // to convert NASA values from days to years. 

    SolarSystem solarSystem;
    
    // NASA add masses too.
    // Sun
    vec3 sun_pos = { -5.919155772743485E-04,  7.385542274920093E-03, -6.154229355980155E-05};
    vec3 sun_vel = {-7.757658457453950E-06,  2.141834870050543E-06,  1.962698872807737E-07};
    sun_vel *= time_normalization; // since [v] = [au/day]
    double sun_mass = 1.0; // mass in solar masses

    // Mercury
    vec3 mercury_pos = {-2.830259311829709E-01, 1.999944216677920E-01, 4.158711313644133E-02};
    vec3 mercury_vel = {-2.161187944929931E-02, -2.207233009714748E-02, 1.783325912551766E-04};
    mercury_vel *= time_normalization;
    double mercury_mass = (3.3/2)*1e-7;

//    // Venus
//    vec3 venus_pos = {-1.979804002020823E-01,  6.983650310789559E-01, 2.081035900640404E-02};
//    vec3 venus_vel = {-1.952418772472886E-02, -5.665705836375252E-03,  1.048675733150797E-03};
//    venus_vel *= time_normalization; // since [v] = [au/day]
//    double venus_mass = (4.9/2)*1e-6;
//
//    // Earth
//    vec3 earth_pos = {7.568866818711384E-01, 6.483949789024712E-01, -1.019825122556535E-04};
//    vec3 earth_vel = {-1.140300222458833E-02, 1.307463198392195E-02, -8.646472327809290E-07};
//    earth_vel *= time_normalization; // since [v] = [au/day]
//    double earth_mass = 3e-6;
//
//    // Mars
//    vec3 mars_pos = {1.388769709785988E+00, 1.556244779300845E-01, -3.104556796346191E-02};
//    vec3 mars_vel = {-9.586334259664849E-04, 1.511271999086297E-02, 3.401500099368581E-04};
//    mars_vel *= time_normalization; // since [v] = [au/day]
//    double mars_mass = 3.3e-7;
//
//    // Jupiter
//    vec3 jupiter_pos = {-2.530234099776776E+00, -4.724912151622468E+00, 7.619427020895447E-02};
//    vec3 jupiter_vel = {6.563190638070959E-03, -3.202416923549941E-03, -1.334759412082958E-04};
//    jupiter_vel *= time_normalization; // since [v] = [au/day]
//    double jupiter_mass = (1.9/2)*1e-3;
//
//    // Saturn
//    vec3 saturn_pos = {  1.850242332516125E+00, -9.881522337036836E+00,  9.816416616683857E-02};
//    vec3 saturn_vel = {5.176955755662489E-03,  1.009140386167208E-03, -2.236909273787052E-04};
//    saturn_vel *= time_normalization; // since [v] = [au/day]
//    double saturn_mass = (5.5/2)*1e-4;
//
//    // Uranus
//    vec3 uranus_pos = {1.705839546641438E+01,  1.018354541295799E+01, -1.831717748605908E-01};
//    vec3 uranus_vel = {-2.044898236053100E-03,  3.193771066110533E-03,  3.841503574618114E-05};
//    uranus_vel *= time_normalization; // since [v] = [au/day]
//    double uranus_mass = 4.4e-5;
//
//    // Neptune
//    vec3 neptune_pos = { 2.896552399030275E+01, -7.545409155038740E+00, -5.121563515019911E-01};
//    vec3 neptune_vel = {7.701585398168282E-04,  3.056586272651174E-03, -8.029789856472037E-05};
//    neptune_vel *= time_normalization; // since [v] = [au/day]
//    double neptune_mass = (1.03/2)*1e-4;

    // Use with: solarSystem.createCelestialBody( position, velocity, mass );
    solarSystem.createCelestialBody( vec3(0,0,0), vec3(0,0,0), sun_mass, "sun");// stationary sun
//    solarSystem.createCelestialBody( sun_pos, sun_vel, sun_mass, "sun");// mobile sun
    solarSystem.createCelestialBody( vec3(0.3075, 0, 0), vec3(0, 12.44, 0), mercury_mass, "mercury" ); // perihelion distance
//    solarSystem.createCelestialBody( mercury_pos, mercury_vel, mercury_mass, "mercury");
 //   solarSystem.createCelestialBody( venus_pos, venus_vel, venus_mass, "venus");
 //   solarSystem.createCelestialBody( earth_pos, earth_vel, earth_mass, "earth");
 //   solarSystem.createCelestialBody( mars_pos, mars_vel, mars_mass, "mars");
 //   solarSystem.createCelestialBody( jupiter_pos, jupiter_vel, jupiter_mass, "jupiter");
 //   solarSystem.createCelestialBody( saturn_pos, saturn_vel, saturn_mass, "saturn");
 //   solarSystem.createCelestialBody( uranus_pos, uranus_vel, uranus_mass, "uranus");
 //   solarSystem.createCelestialBody( neptune_pos, neptune_vel, neptune_mass, "neptune");
    solarSystem.calculateForcesAndEnergy();

    // To get a list (a reference, not copy) of all the bodies in the solar system, we use the .bodies() function
    vector<CelestialBody> &bodies = solarSystem.bodies();

    for(unsigned int i = 0; i<bodies.size(); i++) 
    {
        CelestialBody &body = bodies[i]; // Reference to this body
        cout << "The position of "+body.name+" is" << body.position << " with velocity " << body.velocity << endl;
    }
    
    // setup recording devices
    solarSystem.calculateForcesAndEnergy();
    vec3 angularmomentum;
    angularmomentum = solarSystem.angularMomentum();
    //for perihelion
    double rmin0, rmin1;
    vec3 pos_perihelion0, pos_perihelion1;
    rmin0 = 12;
    rmin1 = 12;

    // save or not
    string answer;
    cout << "do you want to printi(y/n)?\n";
    cin >> answer;
    string runName;
    if( answer == "y")
    {
        cout << "name the run: ";
        cin >> runName; 
        cout << "thank you." << endl;
    }
    else
    {
        cout << "calculating." << endl;
    }

    //integrator
    double dt = (double)timelength/(double)numTimesteps;
    Solver integrator(dt);
    for(int timestep=0; timestep<numTimesteps; timestep++) 
    {
        integrator.Verlet(solarSystem);
        //recording perihelion
        if(  bodies[1].position.length() <= rmin0 and (dt*timestep) <= 1 )
        {
            cout << "Boop! \n"; 
            rmin0 = bodies[1].position.length();
            pos_perihelion0 = bodies[1].position;
        }
        if( bodies[1].position.length() <= rmin1 and (dt*timestep) >= 99 )
        {
            cout << "Beep! \n";
            rmin1 = bodies[1].position.length();
            pos_perihelion1 = bodies[1].position;
        }
        if( timestep%(int)1e8 == 0 and answer == "y")
        {
            cout << "timestep: " << timestep << endl;
            solarSystem.writeToFile(runName);
        }
    }
    solarSystem.calculateForcesAndEnergy();
    double epsilon = 0.1;
    double error = (angularmomentum - solarSystem.angularMomentum()).length()/angularmomentum.length();
    if( error  < epsilon)
    {
        cout << "initial angular momentum " << angularmomentum << " and afterwards " << solarSystem.angularMomentum() << "\n";
        cout << "angular momentum is conserved! with a relative change of " << error << " over " << timelength << " years.\n";
    }
    else
    {
        cout << "initial angular momentum " << angularmomentum << " and afterwards " << solarSystem.angularMomentum() << "\n";
        cout << "angular momentum had a relative drift of  " << error << " over " << timelength << " years.\n";
    }


    double thetap0 = atan2(pos_perihelion0[1], pos_perihelion0[0])*(180*3600)/M_PI;
    double thetap1 = atan2(pos_perihelion1[1], pos_perihelion1[0])*(180*3600)/M_PI;
    cout << "With newtonian gravity, there is a precession of " << thetap1-thetap0 << " arcseconds in a century." << endl;

    return 0;
}

