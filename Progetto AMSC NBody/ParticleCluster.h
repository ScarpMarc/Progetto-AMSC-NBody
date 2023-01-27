#pragma once
#include <vector>

#include "Particle.h"
/*
	TODO
	Update clusters at each frame
	Update boundaries at each frame
*/


/// <summary>
/// Represents a cluster of particles, used to optimise the calculation of interactions with other particles/clusters.
/// </summary>
/// <remarks>
/// Effectively, a cluster acts as a virtual particle and retains all its properties. 
/// Methods that make use of Particle's properties (position, speed, acceleration, mass) are the exact same as the base class;
///		each getter method now returns the sum over all particles of the relevant property, which is calculated upon creation or
///		in general when particles are added/removed.
/// Clusters can have either particles OR other clusters among their children.
/// </remarks>
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
		local_max_boundary = Particle<dim>::max_boundary;
		local_min_boundary = Particle<dim>::min_boundary;

		active = false;

		++maxID;
	}

	// TODO
	/// <summary>
	/// Constructor that takes the particle list, computes the sum of each property and initialises the relevant fields.
	/// </summary>
	ParticleCluster(const std::vector < std::shared_ptr<Particle<dim>>>& children) : Particle<dim>(), ID(maxID), children(children)
	{
		for (const std::shared_ptr<Particle<dim>>& p : children)
		{
			Particle<dim>::mass += p->getMass();
			Particle<dim>::pos += p->get_position();
			Particle<dim>::speed += p->get_speed();
			Particle<dim>::accel += p->get_acc();
		}

		local_max_boundary = Particle<dim>::max_boundary;
		local_min_boundary = Particle<dim>::min_boundary;

		active = true;

		++maxID;
	}

	constexpr inline size_t get_children_num() const { return children.size(); }

	constexpr inline size_t get_children_num_recursive() const;

	// Methods that calculate speed, acceleration, update force... and all getters are those of the base class

	// TODO
	//void saveToFile(std::ofstream& outfile) const;
	//TODO
	//void loadFromFile(std::ifstream& infile);
	//TODO
	/// <summary>
	/// Print particle information
	/// </summary>
	//void print() const;

	/// <summary>
	/// Adds a particle to the cluster.
	/// </summary>
	/// <remarks>
	/// Internally, this function checks against <see cref="max_children"/> whether the maximum amount
	///		of children has been exceeded; if so, the cluster will be divided into sub-clusters (usually 4,
	///		see <see cref="max_children"/>
	/// </remarks>
	/// <param name="p">Pointer to the particle to be added</param>
	void add_particle(const std::shared_ptr<Particle<dim>>& p);

	constexpr inline bool isCluster() const { return false; }

	constexpr void update_boundaries();

	inline const Vector<dim>& get_max_boundary() const { return local_max_boundary; }

	inline const Vector<dim>& get_min_boundary() const { return local_min_boundary; }

	/// <summary>
	/// Queries whether the cluster is active
	/// </summary>
	inline bool is_active() { return is_active; }

private:
	static unsigned int maxID;
	unsigned int ID; // cluster id number

	static unsigned int max_children;

	/// <summary>
	/// A cluster is active if it has children and thus must be included in the computation
	/// </summary>
	bool is_active;

	Vector<dim> local_max_boundary;
	Vector<dim> local_min_boundary;

	/// <summary>
	/// Children particles of this object. We assume that all children are either
	///		Particles or ParticleClusters.
	/// </summary>
	std::vector < std::shared_ptr<Particle<dim>>> children;
};

template<unsigned int dim>
constexpr inline size_t ParticleCluster<dim>::get_children_num_recursive() const
{
	size_t output = children.size();
	for (const std::shared_ptr<Particle<dim>>& p : children)
	{
		if (ParticleCluster<dim>* c = dynamic_cast<ParticleCluster<dim>*>(p.get()))
			output += c->get_children_num_recursive();
		else break;
	}
	return output;
}

template<unsigned int dim>
constexpr void ParticleCluster<dim>::update_boundaries()
{
	for (unsigned int d = 0; d < dim; ++d)
	{
		double lmaxb_component = local_max_boundary[d];
		double lminb_component = local_min_boundary[d];
		//double lmaxb_s_comp = Particle<dim>::max_boundary[d];
#pragma omp parallel for shared(lmaxb_component, lminb_component) reduction(max: lmaxb_comp, lminb_component)
		for (unsigned int i = 0; i < children.size(); ++i)
		{
			double this_pos_component = children[i]->get_position()[d];
			lmaxb_component = this_pos_component > lmaxb_component ? this_pos_component : lmaxb_component;
			lminb_component = this_pos_component < lminb_component ? this_pos_component : lminb_component;
		}
		local_max_boundary[d] = lmaxb_component;
		local_min_boundary[d] = lminb_component;
		if (Particle<dim>::max_boundary[d] > lmaxb_component) local_max_boundary[d] = lmaxb_component;
		if (Particle<dim>::min_boundary[d] < lminb_component) local_min_boundary[d] = lminb_component;
	}
}