#include "gravitationalforce.h"

GravitationalForce::GravitationalForce():
{
}

vec3 Newtonian(CelestialBody body1, CelestialBody body2, vec3 R12, double dist)
{
    return -body1.mass*4*M_PI^2*body2.mass*R12/(dist*dist*dist)
}
