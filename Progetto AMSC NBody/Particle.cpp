#include "Particle.h"

template<unsigned int dim>
Vector<dim> Particle<dim>::calcNewPosition(const unsigned int& delta_time)
{
	_updateSpeed(delta_time);
	_updatePos(delta_time);
	return pos;
}

template<unsigned int dim>
void Particle<dim>::updateResultingForce(const Vector<dim>& resulting_force)
{
	for (int i = 0; i < dim; ++i)
	{
		accel[i] = resulting_force[i] / mass;
	}
}

template<unsigned int dim>
void Particle<dim>::_updateSpeed(const unsigned int& delta_time)
{
	for (int i = 0; i < dim; ++i) 
	{
		speed[i] = accel[i] * delta_time;
	}
}

template<unsigned int dim>
void Particle<dim>::_updatePos(const unsigned int& delta_time)
{
	for (int i = 0; i < dim; ++i)
	{
		pos[i] = speed[i] * delta_time;
	}
}
