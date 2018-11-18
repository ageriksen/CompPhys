#include "solarsystem.h"
#include <iostream>
using namespace std;

SolarSystem::SolarSystem() :
    m_kineticEnergy(0),
    m_potentialEnergy(0)
{
}

CelestialBody& SolarSystem::createCelestialBody(vec3 position, vec3 velocity, double mass, string name) {
    m_bodies.push_back( CelestialBody(position, velocity, mass, name) );
    return m_bodies.back(); // Return reference to the newest added celstial body
}

void SolarSystem::calculateForcesAndEnergy()
{
    m_kineticEnergy = 0;
    m_potentialEnergy = 0;
    m_angularMomentum.zeros();
    m_c = 173*365; // speed of light in AU/yr

    for(CelestialBody &body : m_bodies) {
        // Reset forces on all bodies
        body.force.zeros();
    }
    
    for(int i=0; i<numberOfBodies(); i++) 
    {
        vec3 Force;
        CelestialBody &body1 = m_bodies[i];
        for(int j=i+1; j<numberOfBodies(); j++) 
        {
            CelestialBody &body2 = m_bodies[j];
            vec3 deltaRVector = body1.position - body2.position;
            double dr = deltaRVector.length();
            // Calculate the force and potential energy here
            m_angularMomentum += body1.mass*(body1.velocity.cross(body1.position)) + body2.mass*(body2.velocity.cross(body2.position));
            Force = (-4*M_PI*M_PI*body1.mass*body2.mass*deltaRVector/(dr*dr*dr));
                //*(1 + 3*(deltaRVector.cross(body2.velocity)).lengthSquared()/(dr*dr*m_c*m_c));
            //body1.force += Force;
            body2.force -= Force;
            m_potentialEnergy -= 4*pow(M_PI, 2)*body1.mass*body2.mass/dr;
        }

        m_kineticEnergy += 0.5*body1.mass*body1.velocity.lengthSquared();
    }
}

int SolarSystem::numberOfBodies() const
{
    return m_bodies.size();
}

double SolarSystem::totalEnergy() const
{
    return m_kineticEnergy + m_potentialEnergy;
}

double SolarSystem::potentialEnergy() const
{
    return m_potentialEnergy;
}

double SolarSystem::kineticEnergy() const
{
    return m_kineticEnergy;
}

void SolarSystem::writeToFile(string runName)
{
    for( CelestialBody & body : m_bodies)
    {
        std::ofstream file;
        file.open(runName + body.name + ".dat", ofstream::out | ofstream::app);
        file << body.position.x() << " " << body.position.y() << " " << body.position.z() << "\n";
        file.close(); 
    }
}

vec3 SolarSystem::angularMomentum() const
{
    return m_angularMomentum;
}

std::vector<CelestialBody> &SolarSystem::bodies()
{
    return m_bodies;
}
