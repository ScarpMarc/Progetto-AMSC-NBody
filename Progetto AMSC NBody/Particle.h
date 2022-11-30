#pragma once

#include "Types.h"
#include "Vector.h"

/// <summary>
/// Represents the basic particles that interact in the world
/// </summary>
/// <typeparam name="dim">Dimension (1D, 2D, 3D...).</typeparam>
template <unsigned int dim>
class Particle
{
public:

	/// <summary>
	/// Constructor
	/// </summary>
	/// <param name =
	Particle(unsigned int id,
		Vector<dim> position,
		Vector <dim> speed,
		Vector <dim> acceleration,
		double mass) :
		ID(id),
		pos(position),
		speed(speed),
		accel(acceleration),
		mass(mass)
	{}

	/// <summary>
	/// Calculates the force exerted by the other particle on this particle.
	/// </summary>
	/// <param name="other">Other particle</param>
	Vector<dim> calcForce(const Particle& other) const;
	/// <summary>
	/// Calculates the absolute distance between this particle and the other.
	/// </summary>
	/// <param name="other">Other particle</param>
	Vector<dim> calcDistance(const Particle& other) const;
	/// <summary>
	/// Calculates the new particle position based on speed, acceleration and previous position
	/// </summary>
	/// <param name="delta_time">Delta time, expressed in unit time</param>
	/// <returns>The new position</returns>
	Vector<dim> calcNewPosition(const unsigned int& delta_time);

	/// <summary>
	/// Updates this particle's force from the global matrix.
	/// </summary>
	/// <param name="resulting_force">Resulting force from the matrix</param>
	void updateResultingForce(const Vector<dim>& resulting_force);

	/// <summary>
	/// Gets ID number
	/// </summary>
	unsigned int get_particle_id() { return ID; }

	/// <summary>
	/// Gets position, speed and acceleration
	/// </summary>
	Vector<dim> get_position() { return pos; }
	Vector<dim> get_speed() { return speed; }
	Vector<dim> get_acc() { return accel; }

private:
	void _updateSpeed(const unsigned int& delta_time);
	void _updatePos(const unsigned int& delta_time);

	unsigned int ID; // particle id number

	Vector<dim> pos;
	Vector<dim> speed;
	Vector<dim> accel;
	const double mass;

	// TODO Set value
	// Mass constant k
	const double mass_constant_k = 1.0;

	// TODO Electric constant

	// TODO Magnetic constant (?)
};


