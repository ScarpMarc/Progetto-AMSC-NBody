#pragma once
#include <vector>

#include "Particle.h"

/// <summary>
/// Represents a cluster of particles, used to optimise the calculation of interactions with other particles/clusters.
/// Effectively, a cluster acts as a virtual particle and retains all its properties. 
/// Methods that make use of Particle's properties (position, speed, acceleration, mass) are the exact same as the base class;
///		each getter method now returns the sum over all particles of the relevant property, which is calculated upon creation or
///		in general when particles are added/removed.
/// Clusters can have other clusters among their children.
/// </summary>
/// <typeparam name="dim">Dimension of the space this cluster is living in (2D, 3D)</typeparam>
template<unsigned int dim>
class ParticleCluster :
	public Particle<dim>
{
public:
	/// <summary>
	/// Default constructor which initialises all properties to 0, same as the default constructor of Particle.
	/// </summary>
	ParticleCluster() : Particle<dim>(), ID(maxID)
	{
		++maxID;
	}

	//TODO
	/// <summary>
	/// Constructor that takes the particle list, computes the sum of each property and initialises the relevant fields.
	/// </summary>
	/*ParticleCluster(const std::vector < std::shared_ptr<Particle<dim>>> &children) : Particle<dim>(), ID(maxID)
	{
		++maxID;
	}*/

	// Mehtods that calculate speed, acceleration, update force... and all getters are those of the base class

	// TODO
	//void saveToFile(std::ofstream& outfile) const;
	//TODO
	//void loadFromFile(std::ifstream& infile);
	//TODO
	/// <summary>
	/// Print particle information
	/// </summary>
	//void print() const;

	constexpr inline bool isCluster() const { return false; }


private:
	static unsigned int maxID;
	unsigned int ID; // cluster id number

	std::vector < std::shared_ptr<Particle<dim>>> children;
};