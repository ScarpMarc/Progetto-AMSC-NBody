#pragma once

#include <vector>

#include "Particle.h"

template <unsigned int dim>
int do_simulation_step(std::vector<Particle<dim>>& particles, const unsigned int& ticks_step)
{
#pragma omp parallel 
	{
#pragma omp for
		for (int i = 0; i < particles.size(); ++i)
		{
			Vector<dim> temp;
			/*for (unsigned int component = 0; component < dim; ++component)
			{
				double component_force_on_particle = 0; // Automatically initialised to 0
	//#pragma omp parallel for reduction(+: component_force_on_particle)
				for (int j = 1; j < particles.size(); ++j)
				{
					component_force_on_particle += particles[i].calcForce(particles[(i + j) % particles.size()])[component];
				}

				temp[component] = component_force_on_particle;
			}*/
			for (int j = 1; j < particles.size(); ++j)
			{
				temp += particles[i].calcForceCoefficients(particles[(i + j) % particles.size()]);
			}

			particles[i].updateResultingForce(temp * particles[i].getMass() * particles[i].mass_constant_k);
		}
#pragma omp barrier
#pragma omp for
		for (int i = 0; i < particles.size(); ++i)
		{
			particles[i].calcNewPosition(ticks_step);
		}
	}
	return 0;
}