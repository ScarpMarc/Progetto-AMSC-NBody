#pragma once
#include <vector>

#include "Particle.h"
#include <assert.h>
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
	/// Adds a particle to the cluster, or one of its children.
	/// </summary>
	/// <remarks>
	/// Internally, this function checks against whether this cluster already has other clusters as children, 
	///		we will query whether they are active. If so, 
	///		we can send the particle to the appropriate one. If NONE is active, we prune them all and then add 
	///		<c>p</c> directly to the children vector.
	/// Conversely, if we already have other <c>Particle</c>s as children, we will check against
	///		<see cref="max_children_particles"/> whether the maximum amount
	///		has been exceeded: if so, the cluster will be divided into sub-clusters (usually 8,
	///		see <see cref="max_children"/>); if not, we will append <c>p</c> to the children vector.
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
	inline bool is_active() { return active; }

	/// <summary>
	/// Removes all children. Useful when garbage collecting.
	/// </summary>
	//void remove_all_children();

	/// <summary>
	/// Removes all sub-clusters. Useful when garbage collecting.
	/// </summary>
	void remove_all_subclusters();

private:
	static unsigned int maxID;
	unsigned int ID; // cluster id number

	/// <summary>
	/// Creates sub-clusters automatically in order to accomodate <c>n</c> particles.
	/// </summary>
	/// <remarks>
	/// We want to accomodate <c>n</c> particles; this means we may have to crate more nested levels of sub-clusters.
	///		Since this is a costly operation, we will create enough sub-levels to accomodate all particles.
	///	In a perfect world, particles would have a flat distribution within our region and we would get away with
	///		creating <c>log(n)</c> sub-levels; this is not going to happen in reality, and that is why we can call
	///		the garbage collector directly.
	/// In general, this function will be used at the start of the simulation.
	/// </remarks>
	///<param name = "n">Amount of total children to store</param>
	///<param name = "garbage_collect">Whether to remove unused clusters immediately</param>
	void _create_subclusters(const unsigned int& n, const bool& garbage_collect = false);

	/// <summary>
	/// Creates one level of sub-clusters. The amount is specified by <see cref="num_subclusters"/>.
	/// This will assign particles among <c>children</c> to the appropriate sub-cluster.
	/// </summary>
	void _create_subclusters_one_level();

	/// <summary>
	/// Checks if there are active children (assumed to be other <c>ParticleCluster</c>s)
	/// </summary>
	/// <returns>
	/// Whether there are active children
	/// </returns>
	bool _check_has_active_children() const;

	/*void _update_parent_active_children_count()
	{
		if (parent != nullptr)
		{

		}
	}*/

	/// <summary>
	/// Maximum amount of particle children
	/// </summary>
	static unsigned int max_children_particles;

	/// <summary>
	/// Amount of sub-clusters this cluster can be divided into.
	/// MUST be a cube! We will be dividing the space into <c>num_subclusters</c>x<c>num_subclusters</c> sub-clusters.
	/// </summary>
	static unsigned int num_subclusters;

	/// <summary>
	/// A cluster is active if it has children and thus must be included in the computation.
	/// An inactive cluster will be freed up by the garbage collector.
	/// </summary>
	bool active;

	/// <summary>
	/// The amount of children that are active clusters. When a child becomes active or inactive, it informs the parent.
	/// A cluster with only particles as children will have 0 active children.
	/// </summary>
	//unsigned int active_children;

	Vector<dim> local_max_boundary;
	Vector<dim> local_min_boundary;

	/// <summary>
	/// Children particles of this object. We assume that all children are either
	///		Particles or ParticleClusters.
	/// </summary>
	std::vector < std::shared_ptr<Particle<dim>>> children;

	std::shared_ptr<ParticleCluster<dim>> parent;
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
#pragma omp parallel for reduction(max: lmaxb_component, lminb_component)
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

template<unsigned int dim>
void ParticleCluster<dim>::add_particle(const std::shared_ptr<Particle<dim>>& p)
{
	if (children.size() > 0)
	{
		if (!children[0]->isCluster())
		{
			// First child is not a cluster; this means ALL children are particles. 
			if (children.size() == max_children_particles)
			{
				// We are already at maximum capacity -> create children
				_create_subclusters_one_level();
			}
			else
			{
				children.push_back(p);
				active = true;
			}
		}
		else
		{
			// ALL children are other clusters
			if (_check_has_active_children())
			{
				// We have some sub-clusters that are active.

				// TODO
				// Send to appropriate child, mark active if necessary.
				// Sub-clusters will deal with their own children amount
			}
			else
			{
				/* No sub-clusters are active. This means that they store no particles, and they still exist simply because
					the garbage collector has not been called yet.
					We can remove them (and theur children) and add the particle to our own children.
				*/
				remove_all_subclusters();
			}
		}
	}
	else
	{
		children.push_back(p);
	}
}

template<unsigned int dim>
void ParticleCluster<dim>::remove_all_subclusters()
{
	/* We are calling this function either during a garbage collection sweep, or when adding a new
		particle. We have to delete all sub-clusters, and their sub-clusters; no
		particles should caught in the process.
	*/
	for (const std::shared_ptr<Particle<dim>>& p : children)
	{
		assert(p->isCluster());

		dynamic_cast<ParticleCluster<dim>*>(p.get())->remove_all_subclusters();
	}

	children.clear();
}

template<unsigned int dim>
bool ParticleCluster<dim>::_check_has_active_children() const
{
	bool result = false;
	for (const std::shared_ptr<Particle<dim>>& p : children)
	{
		assert(p->isCluster());

		result |= dynamic_cast<ParticleCluster<dim>*>(p.get())->is_active();
	}
	return result;
}
