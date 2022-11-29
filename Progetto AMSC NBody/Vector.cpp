#include "Vector.h"

template <size_t dim>
double Vector<dim>::eu_norm() const
{
	double sum = 0;
	for (unsigned int i = 0; i < dim; ++i) sum += this[i];
	return sqrt(sum);
}