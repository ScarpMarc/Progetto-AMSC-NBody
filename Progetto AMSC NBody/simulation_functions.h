#pragma once

#include <vector>

#include "Particle.h"
#include "ParticleCluster.h"

template <unsigned int dim>
void assign_particle(const Particle<dim>& particle, const ParticleCluster<dim>& root)
{

}

/// <summary>
/// Handles one simulation step.
/// </summary>
/// <remarks> 
/// First, we divide particles into clusters by assigning them to already-existing ones, and then we create new ones if necessary;
///		after that, we query again iteratively until all clusters are (at least) partially full, but not overbooked.
///		We also compute the near-interaction list of each cluster.
/// 
/// After each particle has been assigned to a cluster, we can use the fast-multipole method to compute far-field approximations
///		and compute forces/speeds
/// </remarks> 
/// <typeparam name="dim">Dimension of the problem</typeparam>
/// <param name="particles">Vector containing the particles</param>
/// <param name="ticks_step">Delta ticks time step</param>
/// <returns>0 if no problems occurred</returns>
template <unsigned int dim>
int do_simulation_step(std::vector<Particle<dim>>& particles, const unsigned int& ticks_step)
{
#pragma omp parallel 
	{
#pragma omp for
		for (unsigned int i = 0; i < particles.size(); ++i)
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
		for (unsigned int i = 0; i < particles.size(); ++i)
		{
			particles[i].calcNewPosition(ticks_step);
		}
	}
	return 0;
}