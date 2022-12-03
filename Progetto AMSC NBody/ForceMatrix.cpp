#include "ForceMatrix.h"

template<unsigned int dim>
void ForceMatrix<dim>::removeParticles(const unsigned int& remove_amt)
{
	// TODO
	// !!! Update the actual counters in the end
	// Also check here that it is possible to add. Maybe return int?

	// This must be called to ensure the array subscript vector matches the amount of particles currently in the system.
	_updatePartialSums_remove(remove_amt);
}

template<unsigned int dim>
void ForceMatrix<dim>::addParticles(const unsigned int& add_amt)
{
	// TODO
	// !!! Update the actual counters in the end
	// Also check here that it is possible to add. Maybe return int?

	// This must be called to ensure the array subscript vector matches the amount of particles currently in the system.
	_updatePartialSums_add(add_amt);
}

template<unsigned int dim>
void ForceMatrix<dim>::updateForces(const std::vector<std::unique_ptr<Particle<dim>>>& particleInteractions)
{
	for (unsigned int i = 0; i < particleInteractions.size() - 1; ++i)
	{
		for (unsigned int j = i + 1; j < particleInteractions.size(); ++j)
		{
			_setInteraction(i, j, particleInteractions[i]->calcForce(*particleInteractions[j]));
		}
	}
}

template<unsigned int dim>
Vector<dim> ForceMatrix<dim>::getTotalForceOnParticle(const unsigned int& idx) const
{
	// TODO test
	Vector<dim> sum;
	for (unsigned int i = 0; i < current_particle_amt; ++i) sum += this(idx, i);
	return sum;
}

template<unsigned int dim>
void ForceMatrix<dim>::_updatePartialSums_remove(const unsigned int& remove_amt)
{
	partial_sums.resize(current_particle_amt - remove_amt);
}

template<unsigned int dim>
void ForceMatrix<dim>::_updatePartialSums_add(const unsigned int& add_amt)
{
	unsigned int initial_sum;
	if (partial_sums.empty()) initial_sum = 0;
	else initial_sum = *(partial_sums.cend() - 1);
	unsigned int sum = initial_sum;
	for (unsigned int i = ForceMatrix<dim>::current_particle_amt; i < current_particle_amt + add_amt; ++i)
	{
		sum += i;
		partial_sums.push_back(sum);
	}
}


/*{
	unsigned int init_sum;
	if (partial_sums2.empty()) init_sum = 0;// partial_sums2.push_back(0);
	else init_sum = *(partial_sums2.cend() - 1);
	unsigned int sum = init_sum;
	for (unsigned int i = 0 particle_amt; i < 0 + 3; ++i)
	{
		sum += i;
		partial_sums2.push_back(sum);
	}
}*/

/*

// Funzione di test per i subscript
void funct(int row, int col)
{
	int a;
	if (row < col)
	{
		a = (current_particle_amt - 1) * row - partial_sums2[row] + col - 1;
		std::cout << a << std::endl;
	}
	// Values below the diagonal have the - sign but are otherwise equal to the others
	else
	{
		a = ((current_particle_amt - 1) * col - partial_sums2[col] + row - 1);
		std::cout << "-" << a << std::endl;
	}
}
*/