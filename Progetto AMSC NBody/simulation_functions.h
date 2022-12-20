#pragma once

#include <vector>

#include "Particle.h"

template<unsigned int dim>
int do_simulation_step(std::vector<Particle<dim>>& particles, const unsigned int & ticks_step)
{
	for (unsigned int i = 0; i < particles.size(); ++i)
	{
		Vector<dim> total_force_on_particle; // Automatically initialised to 0
#pragma omp parallel for shared(total_force_on_particle) reduction(+: total_force_on_particle)
		for (unsigned int j = 1; j < particles.size(); ++j)
		{
			total_force_on_particle += particles[i].calcForce(particles[(i+j)%particles.size()]);
		}

		particles[i].updateResultingForce(total_force_on_particle);
	}
#pragma omp parallel
	for (unsigned int i = 0; i < particles.size(); ++i)
	{
		particles[i].calcNewPosition(ticks_step);
	}
	return 0;
}