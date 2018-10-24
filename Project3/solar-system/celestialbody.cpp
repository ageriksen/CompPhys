#include "celestialbody.h"

CelestialBody::CelestialBody(vec3 pos, vec3 vel, double mass_, string objectname) {
    position = pos;
    velocity = vel;
    mass = mass_;
    name = objectname;
}

CelestialBody::CelestialBody(double x, double y, double z, double vx, double vy, double vz, double mass_, string objectname) {
    position = vec3(x,y,z);
    velocity = vec3(vx,vy,vz);
    mass = mass_;
    name = objectname;
}

void CelestialBody::resetForce() {
    force.zeros();
}
