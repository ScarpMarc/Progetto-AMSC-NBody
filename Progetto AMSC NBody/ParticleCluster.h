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
	ParticleCluster() : Particle<dim>(max_cluster_ID), nest_depth(0), grid_pos(std::array<unsigned int, dim> {})
	{
		local_max_boundary = Particle<dim>::max_boundary;
		local_min_boundary = Particle<dim>::min_boundary;

		active = false;

		++max_cluster_ID;
	}

	// TODO
	/// <summary>
	/// Constructor that takes the particle list, computes the sum of each property and initialises the relevant fields.
	/// </summary>
	ParticleCluster(const std::vector < std::shared_ptr<Particle<dim>>>& children) : Particle<dim>(max_cluster_ID),
		children(children), nest_depth(0), grid_pos(std::array<unsigned int, dim> {})
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

		++max_cluster_ID;
	}

	/// <summary>
	/// Constructor that takes the parent as argument, along with the position of the will-be child within its quadrant.
	/// </summary>
	ParticleCluster(const ParticleCluster<dim>& parent, const unsigned int& depth,
		const std::array<unsigned int, dim> grid_pos) : Particle<dim>(max_cluster_ID), nest_depth(depth), grid_pos(grid_pos)
	{
		this->parent = std::make_shared<ParticleCluster>(parent);

		for (unsigned int i = 0; i < dim; ++i)
		{
			double dim_step = (parent.local_max_boundary[i] - parent.local_min_boundary[i]) / (double)num_subclusters_per_dim;
			local_max_boundary[i] = parent.get_min_boundary()[i] + (double)(grid_pos[i] + 1) * dim_step;
			local_min_boundary[i] = parent.get_min_boundary()[i] + (double)(grid_pos[i]) * dim_step;
		}

		active = false;

		++max_cluster_ID;
	}

	~ParticleCluster()
	{
		parent = nullptr;
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

	constexpr inline bool isCluster() const { return true; }

	constexpr void update_boundaries();

	inline const Vector<dim>& get_max_boundary() const { return local_max_boundary; }

	inline const Vector<dim>& get_min_boundary() const { return local_min_boundary; }

	/// <summary>
	/// Queries whether the cluster is active
	/// </summary>
	inline bool is_active() const { return active; }

	/// <summary>
	/// Removes all children. Useful when garbage collecting.
	/// </summary>
	//void remove_all_children();

	/// <summary>
	/// Removes all sub-clusters. Useful when garbage collecting.
	/// </summary>
	void remove_all_subclusters();

	void set_parent(const ParticleCluster& p) { parent = std::make_shared<ParticleCluster<dim>>(p); }

	/// <summary>
	/// The "depth" of this cluster, i.e. how many level this is below the "main" cluster
	///		that encompasses the entire field; the main cluster has depth 0.
	/// </summary>
	unsigned int get_nest_depth() const { return nest_depth; }
	/// <summary>
	/// Grid position of this cluster among the children of its parent, in each component,
	///		counted starting from the minimum.
	/// </summary>
	std::array<unsigned int, dim> get_grid_pos() const { return grid_pos; }

	void print_recursive() const;

private:
	static unsigned int max_cluster_ID;

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
	/// Clusters are ordered by dimension: first, all Xs from small to large, then Ys and Zs.
	///		The way they are accessed is just like a n-D array.
	/// </remarks>
	///<param name = "n">Amount of total children to store</param>
	///<param name = "garbage_collect">Whether to remove unused clusters immediately</param>
	void _create_subclusters(const unsigned int& n, const bool& garbage_collect = false);

	/// <summary>
	/// Creates one level of sub-clusters. The amount is specified by <see cref="num_subclusters"/>.
	/// This will assign particles among <c>children</c> to the appropriate sub-cluster.
	/// </summary>
	/// <remarks>
	/// Clusters are ordered by dimension: first, all Xs from small to large, then Ys and Zs.
	///		The way they are accessed is just like a n-D array.
	/// All possible sub-clusters will be created: the rationale is that this function will be 
	///		called when there are too many particles and we need to subdivide; hence we create all
	///		of the sub-clusters and fill them immediately.
	/// </remarks>
	void _create_subclusters_one_level();

	/// <summary>
	/// Checks if there are active children (assumed to be other <c>ParticleCluster</c>s)
	/// </summary>
	/// <returns>
	/// Whether there are active children
	/// </returns>
	bool _check_has_active_children() const;

	/// <summary>
	/// Locates the subcluster that would take the particle, based on its position.
	/// The cluster is not guaranteed to actually exist.
	/// </summary>
	/// <returns>
	/// The subscript to the cluster in <p>children</p>
	/// </returns>
	unsigned int _locate_subcluster_for_particle(const std::shared_ptr<Particle<dim>>& p) const;

	/// <summary>
	/// The "depth" of this cluster, i.e. how many level this is below the "main" cluster
	///		that encompasses the entire field; the main cluster has depth 0.
	/// </summary>
	unsigned int nest_depth;
	/// <summary>
	/// Grid position of this cluster among the children of its parent, in each component,
	///		counted starting from the minimum.
	/// </summary>
	std::array<unsigned int, dim> grid_pos;

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
	/// Amount of sub-clusters per dimension this cluster can be divided into.
	/// This simply means we will be dividing the space into <c>num_subclusters_per_dim</c>^<c>dim</c> sub-clusters.
	/// </summary>
	static unsigned int num_subclusters_per_dim;

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
	/// Clusters are ordered by dimension: first, all Xs from small to large, then Ys and Zs.
	///		The way they are accessed is just like a n-D array.
	/// Important assumption: if a cluster has sub-clusters, the vector is always populated with the appropriate
	///		amount of pointers; they may be <c>nullpointer</c>s but they must exist.
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
				unsigned int subscript = _locate_subcluster_for_particle(p);
				if (children[subscript] != nullptr)
				{
					dynamic_cast<ParticleCluster<dim>*>(children[subscript].get())->add_particle(p);
				}
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
				unsigned int subscript = _locate_subcluster_for_particle(p);
				/* Check that the appropriate child exists. If it does, send the particle down;
				*	if it does not, create it.
				*/
				if (children[subscript] != nullptr)
				{
					dynamic_cast<ParticleCluster<dim>*>(children[subscript].get())->add_particle(p);
				}
			}
			else
			{
				/* No sub-clusters are active. This means that they store no particles, and they still exist simply because
					the garbage collector has not been called yet.
					We can remove them (and their children) and add the particle to our own children.
				*/
				remove_all_subclusters();
				children.push_back(p);
				//p->set_parent(this);
				active = true;
			}
		}
	}
	else
	{
		children.push_back(p);
		active = true;
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
	for (const std::shared_ptr<Particle<dim>>& p : children)
	{
		assert(p->isCluster());

		if (dynamic_cast<ParticleCluster<dim>*>(p.get())->is_active()) return true;
	}
	return false;
}

template<unsigned int dim>
void ParticleCluster<dim>::_create_subclusters_one_level()
{
	/*
		Assumption: either we have children and they are all particles, or we have no children.
		This means that the <c>children</c> vector is either filled with at least one particle,
			or empty.
	*/
#ifndef NDEBUG
	for (const std::shared_ptr<Particle<dim>>& p : children)
	{
		assert(!p->isCluster());
	}
#endif
	/* If vector is not empty, we must preserve the <c>Particle</c>s in it.
	*	We can copy it regardless, and figure out later if it was empty or not
	*/
	std::vector < std::shared_ptr<Particle<dim>>> old = children;
	children.clear();

	// Stores modulos of each dim (needed for generation later)
	// dim = 2 --> {2, 4}; dim = 3 --> {2, 4, 8}, etc.
	std::array<unsigned int, dim> dim_mods;

	unsigned int mul = 1;
	for (unsigned int i = 0; i < dim; ++i)
	{
		dim_mods[i] = mul;
		mul *= num_subclusters_per_dim;
	}

	for (unsigned int i = 0; i < mul; ++i)
	{
		std::array<unsigned int, dim> this_iteration_ids;
		for (unsigned int j = 0; j < dim; ++j)
		{
			this_iteration_ids[j] = (i / dim_mods[j]) % num_subclusters_per_dim;
		}
		children.push_back(std::make_shared<ParticleCluster<dim>>(*this, nest_depth + 1, this_iteration_ids));
	}

	// TODO optimise
	/*
		Add old particles to the new appropriate child
	*/
	for (const std::shared_ptr<Particle<dim>>& p : old)
	{
		unsigned int subscript = _locate_subcluster_for_particle(p);
		dynamic_cast<ParticleCluster<dim>*>(children[subscript].get())->add_particle(p);
	}
}

template<unsigned int dim>
unsigned int ParticleCluster<dim>::_locate_subcluster_for_particle(const std::shared_ptr<Particle<dim>>& p) const
{
	/*
	Check my local boundaries and assign the particle to the appropriate region
*/
#ifndef NDEBUG
	for (unsigned int i = 0; i < dim; ++i)
	{
		assert(p->get_position()[i] >= local_min_boundary[i] && p->get_position()[i] <= local_max_boundary[i]);
	}
#endif
	std::array<unsigned int, dim> candidate_grid_pos;
	Vector cluster_size = local_max_boundary - local_min_boundary;
	for (unsigned int i = 0; i < dim; ++i)
	{
		double step0 = ((p->get_position()[i] - local_min_boundary[i]) / cluster_size[i]) // Get number between 0 and 1
			* num_subclusters_per_dim;
		candidate_grid_pos[i] = (unsigned int)(step0); // Get a number and truncate it
	}

#ifndef NDEBUG
	for (unsigned int i = 0; i < dim; ++i)
	{
		assert(candidate_grid_pos[i] >= 0 && candidate_grid_pos[i] < num_subclusters_per_dim);
	}
#endif
	// Crude power function
	unsigned int subscript = 0;
	for (unsigned int i = 0; i < dim; ++i)
	{
		unsigned int mul = 1;
		for (unsigned int j = dim - i; j < dim; ++j)
		{
			mul *= num_subclusters_per_dim;
		}
		subscript += mul * candidate_grid_pos[i];
	}

#ifndef NDEBUG
	unsigned int _mul = 1;
	for (unsigned int j = 0; j < dim; ++j)
	{
		_mul *= num_subclusters_per_dim;
	}

	for (unsigned int i = 0; i < dim; ++i)
	{
		assert(subscript < _mul&& subscript < children.size());
	}
#endif
	return subscript;
}


template<unsigned int dim>
void ParticleCluster<dim>::print_recursive() const
{
	for (unsigned int i = 0; i < nest_depth; ++i) std::cout << "\t";
	std::cout << "Cluster - ID " << this->get_particle_id() << ", " << (active ? "ACTIVE" : "INACTIVE") << ", ";

	if (children.size() > 0)
	{
		std::cout << children.size() << (children.size() > 1 ? " children." : " child.") << std::endl;
		for (const std::shared_ptr<Particle<dim>>& p : children)
		{
			if (p->isCluster())
			{
				dynamic_cast<ParticleCluster<dim>*>(p.get())->print_recursive();
			}
			else
			{
				for (unsigned int i = 0; i < nest_depth + 1; ++i) std::cout << "\t";
				std::cout << "Particle ID " << p->get_particle_id() << std::endl;
			}
		}
	}
	else
	{
		for (unsigned int i = 0; i < nest_depth; ++i) std::cout << "\t";
		std::cout << "no children." << std::endl;
	}
}
