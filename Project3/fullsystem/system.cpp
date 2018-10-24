#include "system.h"
#include <iostream>

System::System()
{
}

//CelestialBody & System::createCelestialBody(
//        vec3 position, vec3 velocity, double mass, string name)
//{
//    bodies.push_back( CelestialBody(position, velocity, mass, name ) );
//    return bodies.back(); // reference to newest celestial body
//}

void System::resetForces() 
{
    /*
     * Ensures the force vector in each object is set to zero.
     */
    for (CelestialBody *obj : bodies) 
    {
        obj->force = {0,0,0};
    }
}

vec3 SolarSystem::angularMomentum() const
{
    return m_angularMomentum;
}

std::vector<CelestialBody> & System::bodies()
{
    return bodies;
}
