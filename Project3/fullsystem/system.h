#ifndef SYSTEM_H
#define SYSTEM_H

#include <vector>
#include "celestialbody.h"

using std::vector;

class System
{
public: // keep everything in puplic for now.

    // Default constructor
    System();
    // Destructor
    ~System() {} // {} makes explicit statement unnecessary in .cpp file.

//    CelestialBody & createCelestialBody( vec3 position, vec3 velocity, double mass, string name );

    void addObject(CelestialBody* newBody) { bodies.push_back(newBody); }

    void resetForces();


    // Vector of pointers i.e. memory adresses of solar system objects of type CelestialBody
    std::vector<CelestialBody*> bodies;
    double m_kineticEnergy=0;
    double m_potentialenergy=0;
};

#endif // SYSTEM_H
