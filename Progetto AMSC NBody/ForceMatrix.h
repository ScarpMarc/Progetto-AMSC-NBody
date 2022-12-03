#pragma once
#include <vector>
#include <set>

#include "Utils.h"

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
	ForceMatrix(const unsigned int& starting_dim) : current_particle_amt(starting_dim), force_matrix(starting_dim) 
	{
		_updatePartialSums_add(starting_dim); // In the future, call addParticles(starting_dim)
	}

	// TODO
	void addParticles(const unsigned int & add_amt);

	void removeParticles(const unsigned int& add_amt);

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
	Vector<dim> getTotalForceOnParticle(const unsigned int& idx) const;

	/// <summary>
	/// Returns a const reference to the desired matrix element.
	/// </summary>
	/// <param name="row">Row index, starting from 0</param>
	/// <param name="col">Column index, starting from 0</param>
	/// <returns>The value referenced by row and col</returns>
	const Vector<dim>& operator()(const unsigned int& row, const unsigned int& col) const
	{
#ifdef _DEBUG
		if (row >= current_particle_amt || col >= current_particle_amt) throw "ERROR! Requested rows and/or cols out of range.";
#endif

		// A particle exerts no force on itself
		if (row == col) return Vector<dim>();
		// Values above the diagonal have the + sign
		if (row < col) return force_matrix[(current_particle_amt - 1) * row - partial_sums[row] + col - 1];
		// Values below the diagonal have the - sign but are otherwise equal to the others
		else return -force_matrix[(current_particle_amt - 1) * col - partial_sums[col] + row - 1];
	}


private:
	inline void _setInteraction(const unsigned int& row, const unsigned int& col, const Vector<dim>& force)
	{
#ifdef _DEBUG
		if (row >= current_particle_amt || col >= current_particle_amt) throw "ERROR! Requested rows and/or cols out of range.";
#endif
		// If row = col we simply do nothing since a particle cannot exert force on itself
		if (row < col) force_matrix[(current_particle_amt - 1) * row - partial_sums[row] + col - 1] = force;
		else if (row < col) force_matrix[(current_particle_amt - 1) * col - partial_sums[col] + row - 1] = force;
	}

	unsigned int current_particle_amt;
	//std::vector<std::unique_ptr<Particle<dim>>> activeParticles;
	// We only keep positive entries of the symmetric matrix (i.e. those above the diagonal)
	std::vector<Vector<dim>> force_matrix;

	/// <summary>
	/// Updates the partial sum vector for quick array subscript computation. MUST be called when adding particles.
	/// </summary>
	/// <param name="add_amt">Amount of particles added</param>
	void _updatePartialSums_add(const unsigned int& add_amt);

	/// <summary>
	/// Updates the partial sum vector for quick array subscript computation. MUST be called when adding particles.
	/// </summary>
	/// <param name="add_amt">Amount of particles added</param>
	void _updatePartialSums_remove(const unsigned int& remove_amt);

	/// <summary>
	/// Keeps a set of pre-computed values to speed up subscript generation in <see cref="operator()"/>.
	/// </summary>
	std::vector<unsigned int> partial_sums;
};