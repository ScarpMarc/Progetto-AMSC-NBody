#include "ForceMatrix.h"

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

// Funzione di test per i subscript
/*void funct(int row, int col)
{
	int a;
	if (row < col)
	{
		a = (current_particle_amt - 1) * row - sumUpTo(1, row+1) + col - 1;
		std::cout << a << std::endl;
	}
	// Values below the diagonal have the - sign but are otherwise equal to the others
	else
	{
		a = ();
		std::cout << "-" << a << std::endl;
	}
}*/