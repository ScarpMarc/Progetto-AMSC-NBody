#include "ForceMatrix.h"

template<unsigned int dim>
void ForceMatrix<dim>::updateForces(const std::vector<std::unique_ptr<Particle<dim>>>& particleInteractions)
{
}

template<unsigned int dim>
Vector<dim> ForceMatrix<dim>::getTotalForceOnParticle(const unsigned int& idx) const
{
	// TODO test
	Vector<dim> sum;
	for (unsigned int i = 0; i < current_particle_amt; ++i) sum += getElementAt(idx, i);
	return sum;
}

template<unsigned int dim>
inline Vector<dim> ForceMatrix<dim>::getElementAt(const unsigned int& row, const unsigned int& col) const
{
	// A particle exerts no force on itself
	if (row == col) return Vector<dim>();
	if (row < col) force_matrix[row * current_particle_amt + col];
	else return -force_matrix[row * current_particle_amt + col];
}
