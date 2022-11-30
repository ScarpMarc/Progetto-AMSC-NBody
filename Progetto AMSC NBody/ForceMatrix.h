#pragma once
#include <vector>
#include <set>

#include "Particle.h"
#include "Types.h"
#include "Vector.h"

/// <summary>
/// Represents the force interactions between all the particles and provides methods for calculating resulting forces that hide parallelism.
/// </summary>
/// <typeparam name="vecDim">Dimension of the vectors</typeparam>
template <unsigned int dim>
class ForceMatrix
{
public:
	ForceMatrix(const unsigned int& starting_dim) : current_particle_amt(starting_dim), force_matrix(starting_dim) {}

	// TODO
	void addParticle();

	// TODO capire qual è il modo più efficiente
	/// <summary>
	/// Updates all forces in the internal matrix. Internally, it queries particles and asks them to calculate their attraction/repulsion to all other particles.
	/// </summary>
	/// <param name="particleInteractions">List of particles</param>
	void updateForces(const std::vector<std::unique_ptr<Particle<dim>>>& particleInteractions);

	/// <summary>
	/// Gets the sum of forces exerted on <param>p</param>
	/// </summary>
	/// <param name="p">The particle of interest</param>
	Vector<dim> getTotalForceOnParticle(const unsigned int &idx) const;

	inline Vector<dim> getElementAt(const unsigned int& row, const unsigned int& col) const;

private:
	unsigned int current_particle_amt;
	//std::vector<std::unique_ptr<Particle<dim>>> activeParticles;
	// We only keep positive entries of the symmetric matrix (i.e. those above the diagonal)
	std::vector<Vector<dim>> force_matrix;
};