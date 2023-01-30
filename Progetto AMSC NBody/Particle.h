#pragma once

#include "Vector.h"
#include "Constants.h"

#include <iomanip>
#include <cmath>
#include <fstream>

/// <summary>
/// Represents the basic particles that interact in the world
/// </summary>
/// <typeparam name="dim">Dimension (1D, 2D, 3D...).</typeparam>
template <unsigned int dim>
class Particle
{
public:
	Particle(
		const unsigned int& ID,
		Vector<dim> position,
		Vector<dim> speed,
		Vector<dim> acceleration,
		double mass) :
		ID(ID),
		pos(position),
		speed(speed),
		accel(acceleration),
		mass(mass)
	{
		for (unsigned int d = 0; d < dim; ++d)
		{
			if (max_boundary[d] < pos[d]) max_boundary[d] = pos[d];
			if (min_boundary[d] > pos[d]) min_boundary[d] = pos[d];
		}

		++maxID;
	}
	/// <summary>
	/// Constructor
	/// </summary>
	/// <param name =
	Particle(
		Vector<dim> position,
		Vector<dim> speed,
		Vector<dim> acceleration,
		double mass) : Particle(maxID, position, speed, acceleration, mass) {}

	Particle(unsigned int ID) : Particle(ID, {}, {}, {}, 0.0) {}

	Particle() : Particle({}, {}, {}, 0.0) {}
	

	/// <summary>
	/// Needed to be able to cast to ParticleCluster through smart pointers.
	/// </summary>
	virtual ~Particle() = default;

	/// <summary>
	/// DEPRECATED Calculates the force exerted by the other particle on this particle.
	/// </summary>
	/// <param name="other">Other particle</param>
	Vector<dim> calcForce(const Particle<dim>& other) const;

	/// <summary>
	/// Calculates the force exerted by the other particle on this particle but omits 
	///		multiplication by the mass constant and the particle mass in order to 
	///		make the calculation more efficient.
	/// </summary>
	/// <typeparam name="dim"></typeparam>
	Vector<dim> calcForceCoefficients(const Particle<dim>& other) const;

	/// <summary>
	/// Calculates the absolute distance between this particle and the other.
	/// </summary>
	/// <param name="other">Other particle</param>
	Vector<dim> calcDistance(const Particle<dim>& other) const;

	/// <summary>
	/// Calculates the new particle position based on speed, acceleration and previous position
	/// </summary>
	/// <param name="delta_time">Delta time, expressed in unit time</param>
	/// <returns>The new position</returns>
	const Vector<dim>& calcNewPosition(const unsigned int& delta_ticks);

	/// <summary>
	/// Updates this particle's force from the global matrix.
	/// </summary>
	/// <param name="resulting_force">Resulting force from the matrix</param>
	void updateResultingForce(const Vector<dim>& resulting_force);

	void saveToFile(std::ofstream& outfile) const;

	void loadFromFile(std::ifstream& infile);

	/// <summary>
	/// Print particle information
	/// </summary>
	void print() const;

	/// <summary>
	/// Gets particle mass
	/// </summary>
	inline const double& getMass() const { return mass; }

	/// <summary>
	/// Gets particle position
	/// </summary>
	inline const Vector<dim>& get_position() const { return pos; }

	/// <summary>
	/// Gets particle speed
	/// </summary>
	inline const Vector<dim>& get_speed() const { return speed; }

	/// <summary>
	/// Gets particle acceleration
	/// </summary>
	inline const Vector<dim>& get_acc() const { return accel; }

	/// <summary>
	/// Gets particle ID
	/// </summary>
	inline const unsigned int& get_particle_id() const { return ID; }

	inline static const Vector<dim>& get_global_max_boundary() { return max_boundary; }

	inline static const Vector<dim>& get_global_min_boundary() { return min_boundary; }

	virtual inline bool isCluster() const { return false; }

	// TODO Set value
	// Mass constant k
	const double mass_constant_k = 6.673e-11;//0.001;

	//void set_parent(const Particle& p) { parent = std::make_shared<Particle>(p); }

protected:
	void _updateSpeed(const unsigned int& delta_ticks);
	void _updatePos(const unsigned int& delta_ticks);

	unsigned int ID; // particle id number

	Vector<dim> pos;
	Vector<dim> speed;
	Vector<dim> accel;
	double mass;

	const double tol = 1e-3;

	static Vector<dim> max_boundary; // Updated by clusters
	static Vector<dim> min_boundary; // Updated by clusters

private:
	static unsigned int maxID;

	// HACK
	//std::shared_ptr<Particle<dim>> parent;

	// TODO Electric constant

	// TODO Magnetic constant (?)
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MEMBER FUNCTIONS
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<unsigned int dim>
void Particle<dim>::print() const
{
	constexpr char coord[3] = { 'X', 'Y', 'Z' };
	std::cout << "Particle ID: " << ID << std::endl;
	std::cout << "Position" << std::endl;
	for (int i = 0; i < dim; ++i)
	{
		std::cout << std::ios::right << std::setw(15) << coord[i] << ": " << pos[i] << std::endl;
	}
	std::cout << "Speed" << std::endl;
	for (int i = 0; i < dim; ++i)
	{
		std::cout << std::ios::right << std::setw(15) << coord[i] << ": " << speed[i] << std::endl;
	}
	std::cout << "Acceleration" << std::endl;
	for (int i = 0; i < dim; ++i)
	{
		std::cout << std::ios::right << std::setw(15) << coord[i] << ": " << speed[i] << std::endl;
	}
	std::cout << std::endl;
}

template<unsigned int dim>
const Vector<dim>& Particle<dim>::calcNewPosition(const unsigned int& delta_ticks)
{
	_updateSpeed(delta_ticks);
	_updatePos(delta_ticks);
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
void Particle<dim>::_updateSpeed(const unsigned int& delta_ticks)
{
	for (int i = 0; i < dim; ++i)
	{
		speed[i] += accel[i] * ((double)delta_ticks / ticks_per_second);
	}
}

template<unsigned int dim>
void Particle<dim>::_updatePos(const unsigned int& delta_ticks)
{
	for (int i = 0; i < dim; ++i)
	{
		pos[i] += speed[i] * ((double)delta_ticks / ticks_per_second);

#pragma omp critical
		{
			if (max_boundary[i] < pos[i]) max_boundary[i] = pos[i];
			if (min_boundary[i] > pos[i]) min_boundary[i] = pos[i];
		}
	}
}

template<unsigned int dim>
Vector<dim> Particle<dim>::calcForceCoefficients(const Particle<dim>& other) const
{
	Vector<dim> displacement = calcDistance(other);
	double distance = displacement.euNorm();
	if (distance <= tol)
	{
		return Vector<dim>();
	}
	return displacement * other.getMass() / (pow(distance, 3));
}

template<unsigned int dim>
Vector<dim> Particle<dim>::calcForce(const Particle<dim>& other) const
{
	return calcForceCoefficients(other) * mass_constant_k * mass;
}

template<unsigned int dim>
Vector<dim> Particle<dim>::calcDistance(const Particle<dim>& other) const
{
	Vector<dim> displacement = Vector<dim>();
	for (int i = 0; i < dim; ++i)
	{
		displacement[i] = other.pos[i] - pos[i];
	}
	return displacement;
}

template<unsigned int dim>
void Particle<dim>::saveToFile(std::ofstream& outfile) const
{
	outfile << pos << speed << accel << mass;
}

template<unsigned int dim>
void Particle<dim>::loadFromFile(std::ifstream& infile)
{
	infile >> pos >> speed >> accel >> mass;
}
