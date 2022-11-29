#pragma once
#include <vector>

#include "Particle.h"
#include "Types.h"

/// <summary>
/// Represents the force interactions between all the particles and provides methods for calculating resulting forces that hide parallelism.
/// </summary>
/// <typeparam name="vecDim">Dimension of the vectors</typeparam>
template <size_t dim>
class ForceMatrix
{
public:
	/// <summary>
	/// Updates all forces in the internal matrix. Internally, it queries particles and asks them to calculate their attraction/repulsion to all other particles.
	/// </summary>
	/// <param name="particleInteractions">List of particles</param>
	void updateForces(std::vector<Particle<dim>> particleInteractions);
	/// <summary>
	/// Gets the sum of forces exerted on <param>p</param>
	/// </summary>
	/// <param name="p">The particle of interest</param>
	Vec<dim> getForceOnParticle(const Particle<dim>& p) const;

private:

};