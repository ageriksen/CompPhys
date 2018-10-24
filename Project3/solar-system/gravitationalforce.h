#ifndef GRAVITATIONALFORCE_H
#define GRAVITATIONALFORCE_H

#include "vec3.h"
#include "celestialbody.h"

using namespace std;

class GravitationalForce
{
public:
    CelestialBody & body1;
    CelestialBody & body2;
    vec3 Newtonian( CelestialBody body1, CelestialBody body2, vec3 R12, double dist);
};

#endif // GRAVITATIONALFORCE_H
