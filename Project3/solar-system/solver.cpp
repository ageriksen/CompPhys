#include "solver.h"
#include "solarsystem.h"

Solver::Solver(double dt) :
    m_dt(dt)
{
}

void Solver::Euler(SolarSystem &system)
{
    system.calculateForcesAndEnergy();

    for(CelestialBody &body : system.bodies()) {
        body.position += body.velocity*m_dt;
        body.velocity += body.force / body.mass * m_dt;
    }
}

void Solver::Verlet(SolarSystem & system)
{
	system.calculateForcesAndEnergy();
	vec3 ai, aii;
	for ( CelestialBody & body: system.bodies() )
	{
		ai = body.force/body.mass;
		body.position += body.velocity*m_dt + 0.5*m_dt*m_dt*ai;
        system.calculateForcesAndEnergy();
		aii = body.force/body.mass;
		body.velocity += 0.5*m_dt*(aii + ai);
	}
}
