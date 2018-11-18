#include "solver.h"
#include "solarsystem.h"

Solver::Solver(double dt) :
    m_dt(dt)
{
}

void Solver::Euler(SolarSystem &system)
{
    /*
     * Forward euler uses the Taylor approx of x and v to find
     * v(t+dt) += a(x(t))*dt
     * x(t+dt) += v(t+dt)*dt
     */
    system.calculateForcesAndEnergy();

    for(CelestialBody &body : system.bodies()) {
        body.velocity += body.force / body.mass * m_dt;
        body.position += body.velocity*m_dt;
    }
}

void Solver::Verlet(SolarSystem & system)
{
    /*
     * Velocity Verlet, making use of the difference
     * between Taylor expansions of v and x to 
     * find next step as 
     * x(t+dt) += v(t)*dt + (1/2)dt^2a(x(t))
     * v(t+1) += (1/2)*dt( a( x(t+dt) ) + a( x(t) ) )
     */
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
