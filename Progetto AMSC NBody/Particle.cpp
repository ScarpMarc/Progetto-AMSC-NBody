#pragma once
#include "Particle.h"

template<size_t dim>
Vector<dim> Particle<dim>::calcNewPosition(const unsigned int& delta_time)
{
	_updateSpeed(delta_time);
	_updatePos(delta_time);
	return pos;
}

template<size_t dim>
inline void Particle<dim>::updateResultingForce(const Vector<dim>& resulting_force)
{
	// TODO dividere ogni elemento del vettore per la massa in parallelo? Si può fare direttamente resulting_force/mass?
}

template<size_t dim>
void Particle<dim>::_updateSpeed(const unsigned int& delta_time)
{
}

template<size_t dim>
void Particle<dim>::_updatePos(const unsigned int& delta_time)
{
}
