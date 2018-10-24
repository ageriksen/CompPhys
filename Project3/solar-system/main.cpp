#include <iostream>
#include <cmath>
#include <cstdlib>
#include "solarsystem.h"
#include "solver.h"
using namespace std;

int main(int numArguments, char **arguments)
{
    int numTimesteps = 2e6;
    int timelength = 100;

    SolarSystem solarSystem;
   
    // Use with: solarSystem.createCelestialBody( position, velocity, mass );
    double earthinitvel = 2*M_PI; // initial velocity circular orbit
    //double earthinitvel = 2*sqrt((double)2)*M_PI; // initial escape velocity
    //double earthinitvel = 1.9*sqrt((double)2)*M_PI; // greater than escape velocity
    double jupiterinitvel = 2*M_PI/5.2;
    solarSystem.createCelestialBody( vec3(0,0,0), vec3(0,0,0), 1.0, "sun");
    //solarSystem.createCelestialBody( vec3(1, 0, 0), vec3(0, earthinitvel, 0), 3e-6, "earth");
    //solarSystem.createCelestialBody( vec3(5.2, 0, 0), vec3(0, jupiterinitvel, 0), (1.9/2)*1e-3, "jupiter");
    //solarSystem.createCelestialBody( vec3(1.52, 0, 0), vec3(0, M_PI/1.52, 0), (1.52/2)*1e-7, "mars" );
    //solarSystem.createCelestialBody( vec3(0.72, 0, 0), vec3(0, M_PI/0.72, 0), (4.9/2)*1e-6, "venus" );
    //solarSystem.createCelestialBody( vec3(9.54, 0, 0), vec3(0, M_PI/9.54, 0), (5.5/2)*1e-4, "saturn" );
    solarSystem.createCelestialBody( vec3(0.3075, 0, 0), vec3(0, 12.44, 0), (3.3/2)*1e-7, "mercury" );
    //solarSystem.createCelestialBody( vec3(19.19, 0, 0), vec3(0, M_PI/19.19, 0), (8.8/2)*1e-5, "uranus" );
    //solarSystem.createCelestialBody( vec3(30.06, 0, 0), vec3(0, M_PI/30.06, 0), (1.03/2)*1e-4, "neptune" );
    //solarSystem.calculateForcesAndEnergy();

    // To get a list (a reference, not copy) of all the bodies in the solar system, we use the .bodies() function
    vector<CelestialBody> &bodies = solarSystem.bodies();
    //bodies[0].velocity = -1*solarSystem.angularMomentum();
    //double Mtot;
    //vec3 rCOM, RCOM; 
    //for( CelestialBody & body: bodies)
    //{
    //    Mtot += body.mass;
    //    rCOM += body.mass*body.position;
    //}
    //RCOM = rCOM/Mtot;
    //bodies[0].position = -1*RCOM;


    for(unsigned int i = 0; i<bodies.size(); i++) 
    {
        CelestialBody &body = bodies[i]; // Reference to this body
        cout << "The position of "+body.name+" is" << body.position << " with velocity " << body.velocity << endl;
    }

    vec3 angularmomentum;
    solarSystem.calculateForcesAndEnergy();
    angularmomentum = solarSystem.angularMomentum();
    double thetap0 = M_PI*atan(bodies[1].position[0]/bodies[1].position[1])/(180*3600);

    string answer;
    cout << "do you want to printi(y/n)?\n";
    cin >> answer;
    double dt = (double)timelength/(double)numTimesteps;
    string runName;
    if( answer == "y")
    {
        cout << "name the run: ";
        cin >> runName; 
    }
    Solver integrator(dt);
    for(int timestep=0; timestep<numTimesteps; timestep++) 
    {
        integrator.Verlet(solarSystem);
        if( timestep%1 == 0 and answer == "y")
        {
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


    double thetap1 = M_PI*(atan(bodies[1].position[0]/bodies[1].position[1]))/(180*3600);

    cout << "With relativistic gravity, there is a precession of " << thetap1-thetap0 << " arcseconds in a century." << endl;

    return 0;
}

