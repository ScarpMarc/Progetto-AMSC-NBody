#pragma once
#include <vector>

#include "Particle.h"
#include "Globals.h"
#include <assert.h>
#include <set>
#include <map>
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
template <unsigned int dim>
class ParticleCluster : public Particle<dim>
{
public:
	/// <summary>
	/// Default constructor which initialises all properties to 0, same as the default constructor of Particle.
	/// </summary>
	ParticleCluster() : Particle<dim>(max_cluster_ID), nest_depth(0), grid_pos(std::array<unsigned int, dim>{}), last_boundary_update(0),
		parent(nullptr), has_particles(false)
	{
		// Brutal hack to avoid having particles exactly on the edge
		local_max_boundary = Particle<dim>::max_boundary + Particle<dim>::max_boundary * .01;
		local_min_boundary = Particle<dim>::min_boundary + Particle<dim>::min_boundary * .01;

		double max = local_max_boundary[0], min = local_min_boundary[0];
		for (unsigned int i = 0; i < dim; ++i)
		{
			if (local_max_boundary[i] > max) max = local_max_boundary[i];
			if (local_min_boundary[i] < min) min = local_min_boundary[i];
		}

		for (unsigned int i = 0; i < dim; ++i)
		{
			local_max_boundary[i] = max;
			local_min_boundary[i] = min;
		}

		active = false;

		++max_cluster_ID;
	}

	// Deprecated! Dangerous!
	/// <summary>
	/// Constructor that takes the particle list, computes the sum of each property and initialises the relevant fields.
	/// </summary>
	/*ParticleCluster(const std::vector<unsigned int>& particle_ids) : Particle<dim>(max_cluster_ID),
		nest_depth(0), grid_pos(std::array<unsigned int, dim>{}), last_boundary_update(0), parent(nullptr)
	{
		for (const unsigned int& i : particle_ids)
		{
			add_particle(i);
		}
		for (const unsigned int& i : children_particles)
		{
			Particle<dim>::mass += _get_particle_global(i).getMass();
			Particle<dim>::pos += _get_particle_global(i).get_position();
			Particle<dim>::speed += _get_particle_global(i).get_speed();
			Particle<dim>::accel += _get_particle_global(i).get_acc();
		}
		local_max_boundary = Particle<dim>::max_boundary;
		local_min_boundary = Particle<dim>::min_boundary;
		active = true;
		++max_cluster_ID;
	}*/

	/// <summary>
	/// Constructor that takes the parent as argument, along with the position of the will-be child within its quadrant.
	/// </summary>
	ParticleCluster(ParticleCluster<dim>* parent, const unsigned int& depth,
		const std::array<unsigned int, dim>& grid_pos) : Particle<dim>(max_cluster_ID), nest_depth(depth), grid_pos(grid_pos),
		last_boundary_update(0), parent(parent), has_particles(false)
	{
		for (unsigned int i = 0; i < dim; ++i)
		{
			double dim_step = (parent->local_max_boundary[i] - parent->local_min_boundary[i]) / (double)num_subclusters_per_dim;
			local_max_boundary[i] = parent->get_min_boundary()[i] + (double)(grid_pos[i] + 1) * dim_step;
			local_min_boundary[i] = parent->get_min_boundary()[i] + (double)(grid_pos[i]) * dim_step;
		}

		active = false;

		++max_cluster_ID;
	}

	constexpr inline size_t get_children_particle_num() const { return children_particles.size(); }

	constexpr inline size_t get_children_clusters_num() const { return children_clusters.size(); }

	size_t get_children_clusters_num_recursive() const;

	size_t get_children_clusters_num_active_recursive() const;

	size_t get_children_particle_num_recursive() const;

	void TEST_update()
	{
		_update_boundaries_recursive();
		_check_particles_recursive();
		_update_physics();
	}

	// Methods that calculate speed, acceleration, update force... and all getters are those of the base class

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
	/// <param name="p">Particle index from the global container</param>
	void add_particle(const unsigned int& p);

	/// <summary>
	/// Attempts to relocate the particle by either storing it or sending it to the parent. Effectively,
	///		an "attempt-to-add" function.
	/// </summary>
	/// <remarks>
	/// When particles are no longer within boundaries, we cannot assume anything about where they went, since
	///		many factors contribute to their movement relative to the cluster. What we can do to find a reasonable
	///		compromise is to pass the particle to the parent, which will check against its boundaries and either
	///		keep it or pass it on (to another child or to its own parent).
	/// This operation may result in some clusters sub-dividing, but at the moment I have no idea how else to implement it.
	/// </remarks>
	/// <param name="p"></param>
	void relocate_particle(const unsigned int& p);

	/// <summary>
	/// Checks all sub-clusters and eliminates those that are no longer active.
	/// </summary>
	void garbage_collect();

	constexpr inline bool isCluster() const { return true; }

	inline const Vector<dim>& get_max_boundary() const { return local_max_boundary; }

	inline const Vector<dim>& get_min_boundary() const { return local_min_boundary; }

	/// <summary>
	/// Queries whether the cluster is active
	/// </summary>
	inline bool is_active() const { return active; }

	/// <summary>
	/// Removes all sub-clusters. Useful when garbage collecting.
	/// </summary>
	void remove_all_subclusters_recursive();

	// void set_parent(ParticleCluster const& p) { parent = p; }

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

	/// <summary>
	/// Crude function to test if the cluster has sub-cluster or not. Will get around to making a better version if I have time;
	///		this function is here to avoid having to refactor the cluster-particle check in the future.
	/// </summary>
	/// <returns>TRUE if the cluster contains particles; FALSE if it contains other clusters.</returns>
	bool contains_particles() const
	{
		return has_particles;
	}

	bool is_empty() const
	{
		return children_clusters.size() == 0 && children_particles.size() == 0;
	}

private:
	static unsigned int max_cluster_ID;

	/// <summary>
	/// Used to update the boundaries, bottom-down from the father.
	/// </summary>
	void _update_boundaries_recursive();

	/// <summary>
	/// Used to check whether all particles are within boundaries, after calling <link cref="update_boundaries_recursive()"/>
	/// If not (they have moved), we relocate them.
	/// </summary>
	void _check_particles_recursive();

	/// <summary>
	/// Calculates center of mass, total speed ecc. in a parallel fashion. Called when we have already set all boundaries and checked the particles.
	/// </summary>
	void _update_physics();

	/// <summary>
	/// Calculates the center of mass of all particles (= <see cref="pos"/> of the cluster)
	/// </summary>
	void _calc_center_of_mass();
	/// <summary>
	/// Calculates the vector sum of all the speeds of the particles (= <see cref="speed"/> of the cluster)
	/// </summary>
	void _calc_speed();
	/// <summary>
	/// Calculates the vector sum of all the accelerations of the particles (= <see cref="accel"/> of the cluster)
	/// </summary>
	void _calc_acc();
	/// <summary>
	/// Calculates sum of all masses of the particles (= <see cref="mass"/> of the cluster)
	/// </summary>
	void _calc_mass();

	/// <summary>
	/// Gathers particles from sub-clusters and appends them to the local vector.
	/// To be used ONLY when relocating. There is no going back!
	/// </summary>
	/// <remarks>
	/// This function takes particles from the sub-clusters, hence it is never called on a cluster
	///		that contains particles itself.
	/// </remarks>
	void __gather_particles_recursive();

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
	///< param name = "n">Amount of total children to store</param>
	///< param name = "garbage_collect">Whether to remove unused clusters immediately</param>
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
	/// Creates a single sub-cluster at a certain region, as specified by <c>location</c>. If there is already
	///		a sub-cluster in that position, the function fails silently.
	/// </summary>
	/// <remarks>
	/// Each index in <c>location</c> must be between 0 and <see cref="num_subclusters_per_dim"/>-1.
	/// </remarks>
	/// <param name = "location">Location of subcluster</param>
	void __create_subcluster_at(const std::array<unsigned int, dim> location);

	/// <summary>
	/// Gets the subscript to the cluster at <c>location</c>. The cluster is not guaranteed to exist.
	/// </summary>
	/// <remarks>
	/// Each index in <c>location</c> must be between 0 and <see cref="num_subclusters_per_dim"/>-1.
	/// </remarks>
	/// <param name = "location">Location of subcluster</param>
	unsigned int _get_cluster_subscript(const std::array<unsigned int, dim> location) const;

	/// <summary>
	/// Locates the subcluster that would take the particle, based on its position.
	/// The cluster is not guaranteed to actually exist.
	/// Effectively equal to calling _get_cluster_subscript(_locate_quadrant_for_particle(p)).
	/// </summary>
	/// <returns>
	/// The subscript to the cluster in <p>children</p>
	/// </returns>
	unsigned int _locate_subcluster_for_particle(const Particle<dim>& p) const;

	/// <summary>
	/// Locates the right quadrant for a particle. That quadrant is not guaranteed to
	///		have a cluster.
	/// </summary>
	std::array<unsigned int, dim> _locate_quadrant_for_particle(const Particle<dim>& p) const;

	/// <summary>
	/// Checks if there are active sub-clusters.
	/// </summary>
	/// <returns>
	/// Whether there are active children
	/// </returns>
	bool _check_has_active_subclusters() const;

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

	inline Particle<dim>& _get_particle_global(const unsigned int& idx) const { return global_particles[idx]; }

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
	// unsigned int active_children;

	Vector<dim> local_max_boundary;
	Vector<dim> local_min_boundary;

	// std::vector<std::shared_ptr<Particle<dim>>> children;
	/// <summary>
	/// Sub-clusters of this cluster. We assume that all children are either
	///		Particles (contained in <see cref="children_particles"/> or ParticleClusters.
	/// </summary>
	/// <remarks>
	/// Clusters are ordered by dimension: first, all Xs from small to large, then Ys and Zs.
	/// The map maps position to cluster.
	/// </remarks>
	std::map<unsigned int, ParticleCluster<dim>> children_clusters;

	/// <summary>
	/// Particles managed by this cluster. We assume that all children are either
	///		Particles ParticleClusters (contained in <see cref="children_clusters"/>.
	/// </summary>
	/// <remarks>
	///
	/// </remarks>
	std::set<unsigned int> children_particles;

	ParticleCluster<dim>* parent;

	bool has_particles;

	unsigned int last_boundary_update;
};

template <unsigned int dim>
size_t ParticleCluster<dim>::get_children_particle_num_recursive() const
{
	size_t output = children_particles.size();
	for (const auto& c : children_clusters)
	{
		output += c.second.get_children_particle_num_recursive();
	}
	return output;
}

template <unsigned int dim>
void ParticleCluster<dim>::relocate_particle(const unsigned int& p)
{
	bool relocated = false;
	const Particle<dim>& ref = _get_particle_global(p);
	for (unsigned int i = 0; i < dim; ++i)
	{
		// If the particle is not within our boundaries
		if (ref.get_position()[i] > local_max_boundary[i] || ref.get_position()[i] < local_min_boundary[i])
		{
			// Relocate
			parent->relocate_particle(p);
			relocated = true;
			break;
		}
	}
	if (!relocated)
	{
		/* If we get here it means that the particle actually was within boundaries.
		 * We can add it to ourselves.
		 */

		this->add_particle(p);
	}
}

template <unsigned int dim>
void ParticleCluster<dim>::_check_particles_recursive()
{
	for (auto& c : children_clusters)
	{
		c.second._check_particles_recursive();
	}

	//#pragma omp single
	for (auto it = children_particles.begin(); it != children_particles.end();)
	{
		bool removed = false;
		for (unsigned int j = 0; j < dim; ++j)
		{
			// If particle
			if (_get_particle_global(*it).get_position()[j] > local_max_boundary[j] || _get_particle_global(*it).get_position()[j] < local_min_boundary[j])
			{
				// relocate
				parent->relocate_particle(*it);
				removed = true;
				break;
			}
		}
		if (!removed)
			++it;
		else
			it = children_particles.erase(it);
	}
}

template <unsigned int dim>
void ParticleCluster<dim>::_update_boundaries_recursive()
{
	/*
		for (unsigned int d = 0; d < dim; ++d)
		{
			double lmaxb_component = local_max_boundary[d];
			double lminb_component = local_min_boundary[d];
			// double lmaxb_s_comp = Particle<dim>::max_boundary[d];
#pragma omp parallel for reduction(max \
								   : lmaxb_component, lminb_component)
			for (unsigned int i = 0; i < children_particles.size(); ++i)
			{
				double this_pos_component = children[i]->get_position()[d];
				lmaxb_component = this_pos_component > lmaxb_component ? this_pos_component : lmaxb_component;
				lminb_component = this_pos_component < lminb_component ? this_pos_component : lminb_component;
			}
			local_max_boundary[d] = lmaxb_component;
			local_min_boundary[d] = lminb_component;
			if (Particle<dim>::max_boundary[d] > lmaxb_component)
				local_max_boundary[d] = lmaxb_component;
			if (Particle<dim>::min_boundary[d] < lminb_component)
				local_min_boundary[d] = lminb_component;
		}
*/
	if (parent == nullptr)
	{
		// We are the main cluster. Update using globals from Particle
		local_max_boundary = Particle<dim>::max_boundary + Particle<dim>::max_boundary * .01;
		local_min_boundary = Particle<dim>::min_boundary + Particle<dim>::min_boundary * .01;

		// Resize to a cube

		double max = local_max_boundary[0], min = local_min_boundary[0];
		for (unsigned int i = 0; i < dim; ++i)
		{
			if (local_max_boundary[i] > max) max = local_max_boundary[i];
			if (local_min_boundary[i] < min) min = local_min_boundary[i];
		}

		for (unsigned int i = 0; i < dim; ++i)
		{
			local_max_boundary[i] = max;
			local_min_boundary[i] = min;
		}
	}
	else
	{
		for (unsigned int i = 0; i < dim; ++i)
		{
			double dim_step = (parent->get_max_boundary()[i] - parent->get_min_boundary()[i]) / (double)num_subclusters_per_dim;
			local_max_boundary[i] = parent->get_min_boundary()[i] + (double)(grid_pos[i] + 1) * dim_step;
			local_min_boundary[i] = parent->get_min_boundary()[i] + (double)(grid_pos[i]) * dim_step;
		}
	}
	// OMP
#pragma omp parallel for
	for (long i = 0; i < children_clusters.size(); ++i)
	{
		auto c = std::next(children_clusters.begin(), i);
		c->second._update_boundaries_recursive();
	}
}

template <unsigned int dim>
void ParticleCluster<dim>::add_particle(const unsigned int& p)
{
	assert(children_particles.find(p) == children_particles.end());
	if (is_empty())
	{
		children_particles.insert(p);
		active = true;
		has_particles = true;
	}
	else
	{
		if (contains_particles())
		{
			assert(!_check_has_active_subclusters());
			/*
				All children are particles.
			*/
			if (children_particles.size() >= max_children_particles)
			{
				// We are already at maximum capacity -> create children
				_create_subclusters_one_level();
				// Clusters created, particles reassigned -> add p
				unsigned int subscript = _locate_subcluster_for_particle(_get_particle_global(p));
				// WARNING we assume that each child exists since we just created them all. Problems may arise when multithreading...
				assert(children_clusters.find(subscript) != children_clusters.end());
				children_clusters[subscript].add_particle(p);
				has_particles = false;
				assert(children_particles.size() == 0);
			}
			else
			{
				children_particles.insert(p);
				active = true;
				has_particles = true;
				assert(children_clusters.size() == 0);
			}
		}
		else
		{
			assert(!contains_particles() && children_particles.size() == 0);
			// All children are other clusters
			if (!_check_has_active_subclusters())
			{
				/*
				If no sub-clusters are active, this means that they store no particles, and they still exist simply because
					the garbage collector has not been called yet.
				We can remove them (and their children) and add the particle to our own children.
				*/
				remove_all_subclusters_recursive();
				children_particles.insert(p);
				active = true;
				has_particles = true;
				assert(children_clusters.size() == 0);
			}
			else
			{
				// ACTIVE SUB-CLUSTERS
				// Assign particle to right sub-cluster
				unsigned int subscript = _locate_subcluster_for_particle(_get_particle_global(p));
				if (children_clusters.find(subscript) == children_clusters.end())
				{
					__create_subcluster_at(_locate_quadrant_for_particle(_get_particle_global(p)));
				}
				children_clusters[subscript].add_particle(p);
				active = true;
				assert(!has_particles);
			}
		}
	}
}

template <unsigned int dim>
void ParticleCluster<dim>::__gather_particles_recursive()
{
	assert(!has_particles);
	//assert(children_clusters.size() != 0);

		// We can use OMP
#pragma omp parallel for
	for (unsigned int c_idx = 0; c_idx < children_clusters.size(); ++c_idx)
	{
		auto c = std::next(children_clusters.begin(), c_idx);

		if (!c->second.contains_particles()) // First gather sub-particles
		{
			c->second.__gather_particles_recursive();
		}

		for (const unsigned int& i : c->second.children_particles)
		{
			assert(children_particles.find(i) == children_particles.end());
			bool relocated = false;
			// If particle is not within our boundaries, forward it to father
			for (unsigned int j = 0; j < dim; ++j)
			{
				if (_get_particle_global(i).get_position()[j] > local_max_boundary[j] || _get_particle_global(i).get_position()[j] < local_min_boundary[j])
				{
					parent->relocate_particle(i);
					relocated = true;
					break;
				}
			}
			// If we did not relocate the particle, it means we can keep it
			if (!relocated)
			{
#pragma omp critical
				children_particles.insert(i);
			}

		}
#pragma omp critical
		c->second.children_particles.clear();
#pragma omp critical
		c->second.active = false;
	}

	children_clusters.clear();
	has_particles = true;
}

template<unsigned int dim>
void ParticleCluster<dim>::remove_all_subclusters_recursive()
{
	/* We are calling this function either during a garbage collection sweep, or when adding a new
		particle. We have to delete all sub-clusters, and their sub-clusters; no
		particle should get caught in the process.
	*/
	assert(get_children_particle_num() == 0);
	// assert(!is_active());
	for (auto& p : children_clusters)
	{
		assert(p.second.get_children_particle_num() == 0);
		assert(!p.second.is_active());

		p.second.remove_all_subclusters_recursive();
	}

	children_clusters.clear();
}

template <unsigned int dim>
bool ParticleCluster<dim>::_check_has_active_subclusters() const
{
	for (const auto& c : children_clusters)
	{
		if (c.second.is_active())
			return true;
	}
	return false;
}

template <unsigned int dim>
void ParticleCluster<dim>::_create_subclusters_one_level()
{
	/*
		Assumption: either we have children and they are all particles, or we have no children.
		This means that the <c>children</c> vector is either filled with at least one particle,
			or empty.
	*/
	assert((contains_particles() && children_clusters.empty()));

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
		__create_subcluster_at(this_iteration_ids);
	}

	/*
		Add old particles to the new appropriate child
	*/
	for (const unsigned int& p : children_particles)
	{
		bool relocated = false;
		// If particle is not within our boundaries, forward it to father
		for (unsigned int j = 0; j < dim; ++j)
		{
			if (_get_particle_global(p).get_position()[j] > local_max_boundary[j] || _get_particle_global(p).get_position()[j] < local_min_boundary[j])
			{
				parent->relocate_particle(p);
				relocated = true;
				break;
			}
		}
		if (!relocated)
		{
			unsigned int subscript = _locate_subcluster_for_particle(_get_particle_global(p));
			children_clusters[subscript].add_particle(p);
		}
	}
	// Delete old ones.
	children_particles.clear();
}

template <unsigned int dim>
void ParticleCluster<dim>::__create_subcluster_at(const std::array<unsigned int, dim> location)
{
#ifndef NDEBUG
	for (unsigned int i = 0; i < dim; ++i)
	{
		assert(location[i] >= 0 && location[i] < num_subclusters_per_dim);
	}

	assert(children_clusters.find(_get_cluster_subscript(location)) == children_clusters.end());

#endif
	// Create new map entry. Some serious dark magic here
	unsigned int newdepth = nest_depth + 1;
	unsigned int subscript = _get_cluster_subscript(location);
	children_clusters.emplace(std::piecewise_construct,
		std::forward_as_tuple(subscript),
		std::forward_as_tuple(this, newdepth, location));
}

/*ParticleCluster(ParticleCluster<dim>* parent, const unsigned int& depth,
	const std::array<unsigned int, dim> grid_pos)*/

template <unsigned int dim>
unsigned int ParticleCluster<dim>::_get_cluster_subscript(const std::array<unsigned int, dim> location) const
{
#ifndef NDEBUG
	for (unsigned int i = 0; i < dim; ++i)
	{
		assert(location[i] >= 0 && location[i] < num_subclusters_per_dim);
	}
#endif
	unsigned int subscript = 0;
	for (unsigned int i = 0; i < dim; ++i) // Crude power function
	{
		unsigned int mul = 1;
		for (unsigned int j = dim - i; j < dim; ++j)
		{
			mul *= num_subclusters_per_dim;
		}
		subscript += mul * location[i];
	}
#ifndef NDEBUG
	unsigned int mul = 1;
	for (unsigned int i = 0; i < dim; ++i)
	{
		mul *= num_subclusters_per_dim;
	}
	assert(subscript < mul);
#endif

	return subscript;
}

template <unsigned int dim>
std::array<unsigned int, dim> ParticleCluster<dim>::_locate_quadrant_for_particle(const Particle<dim>& p) const
{
#ifndef NDEBUG
	for (unsigned int i = 0; i < dim; ++i)
	{
		assert(p.get_position()[i] >= local_min_boundary[i] && p.get_position()[i] <= local_max_boundary[i]);
	}
#endif
	std::array<unsigned int, dim> candidate_grid_pos;
	Vector cluster_size = local_max_boundary - local_min_boundary;
	for (unsigned int i = 0; i < dim; ++i)
	{
		double step0 = ((p.get_position()[i] - local_min_boundary[i]) / cluster_size[i]) // Get number between 0 and 1
			* num_subclusters_per_dim;
		candidate_grid_pos[i] = (unsigned int)(step0); // Get a number and truncate it
	}

#ifndef NDEBUG
	for (unsigned int i = 0; i < dim; ++i)
	{
		assert(candidate_grid_pos[i] >= 0 && candidate_grid_pos[i] < num_subclusters_per_dim);
	}
#endif

	return candidate_grid_pos;
}

template <unsigned int dim>
unsigned int ParticleCluster<dim>::_locate_subcluster_for_particle(const Particle<dim>& p) const
{
	return _get_cluster_subscript(_locate_quadrant_for_particle(p));
}

template <unsigned int dim>
void ParticleCluster<dim>::print_recursive() const
{
	for (unsigned int i = 0; i < nest_depth; ++i)
		std::cout << "\t";
	std::cout << "Cluster - ID " << this->get_particle_id() << ", " << (active ? "ACTIVE" : "INACTIVE") << ", ";
	if (is_empty())
	{
		std::cout << "empty." << std::endl;
	}
	else
	{
		std::cout << "contains ";

		if (contains_particles())
		{
			std::cout << children_particles.size() << (children_particles.size() > 1 ? " particles." : " particle.") << std::endl;
			for (const unsigned int& p : children_particles)
			{
				for (unsigned int i = 0; i < nest_depth + 1; ++i)
					std::cout << "\t";
				std::cout << "Particle ID " << _get_particle_global(p).get_particle_id() << std::endl;
			}
		}
		else
		{
			std::cout << children_clusters.size() << (children_clusters.size() > 1 ? " clusters." : " cluster.") << std::endl;
			for (const auto& p : children_clusters)
			{
				p.second.print_recursive();
			}
		}
	}
}

template <unsigned int dim>
size_t ParticleCluster<dim>::get_children_clusters_num_active_recursive() const
{
	size_t output = 0;
	for (const auto& p : children_clusters)
	{
		if (p.second.is_active())
		{
			++output;
			output += p.second.get_children_clusters_num_active_recursive();
		}
	}
	return output;
}

template <unsigned int dim>
size_t ParticleCluster<dim>::get_children_clusters_num_recursive() const
{
	size_t output = 0;
	for (const auto& p : children_clusters)
	{
		++output;
		output += p.second.get_children_clusters_num_recursive();
	}
	return output;
}

template <unsigned int dim>
void ParticleCluster<dim>::_update_physics()
{
	if (children_particles.size() == 0)
	{
		active = false;
	}
	else
	{
#ifndef _WIN32 // Tasks are unsupported by MSVC
#pragma omp parallel
		{
#pragma omp single
			{
#pragma omp taskgroup
				{
#pragma omp task
#endif
					_calc_speed();
#ifndef _WIN32 // Tasks are unsupported by MSVC
#pragma omp task
#endif
					_calc_acc();
#ifndef _WIN32 // Tasks are unsupported by MSVC
#pragma omp task
#endif
					_calc_mass();
#ifndef _WIN32 // Tasks are unsupported by MSVC
				}
				// Center of mass depends on total mass, we must run it after the others.
#pragma omp task
#endif
				_calc_center_of_mass();
#ifndef _WIN32 // Tasks are unsupported by MSVC
			}
		}
#endif
	}
}

template <unsigned int dim>
void ParticleCluster<dim>::_calc_center_of_mass()
{
	Vector<dim> temp_pos = {};
#ifndef _WIN32 // Custom sum not supported by MSVC
#pragma omp parallel for reduction(VectorSum \
								   : temp_pos)
#endif
	for (unsigned int i = 0; i < children_particles.size(); ++i)
	{
		auto val = std::next(children_particles.begin(), i);
		temp_pos += _get_particle_global(*val).get_position() * _get_particle_global(*val).getMass() / this->mass;
	}
	this->pos = temp_pos;
}

template <unsigned int dim>
void ParticleCluster<dim>::_calc_speed()
{
	Vector<dim> temp_speed = {};
	// Reduction declared in Particle.h
#ifndef _WIN32 // Custom sum not supported by MSVC
#pragma omp parallel for reduction(VectorSum \
								   : temp_speed)
#endif
	for (unsigned int i = 0; i < children_particles.size(); ++i)
	{
		auto val = std::next(children_particles.begin(), i);
		temp_speed += _get_particle_global(*val).get_speed();
	}
	this->speed = temp_speed;
}

template <unsigned int dim>
void ParticleCluster<dim>::_calc_acc()
{
	Vector<dim> temp_accel = {};
	// Reduction declared in Particle.h
#ifndef _WIN32 // Custom sum not supported by MSVC
#pragma omp parallel for reduction(VectorSum \
								   : temp_accel)
#endif
	for (unsigned int i = 0; i < children_particles.size(); ++i)
	{
		auto val = std::next(children_particles.begin(), i);
		temp_accel += _get_particle_global(*val).get_acc();
	}
	this->accel = temp_accel;
}

template <unsigned int dim>
void ParticleCluster<dim>::_calc_mass()
{
	double temp_mass = 0;
	// Reduction declared in Particle.h
#ifndef _WIN32 // Custom sum not supported by MSVC
#pragma omp parallel for reduction(+ \
								   : temp_mass)
#endif
	for (unsigned int i = 0; i < children_particles.size(); ++i)
	{
		auto val = std::next(children_particles.begin(), i);
		temp_mass += _get_particle_global(*val).getMass();
	}
	this->mass = temp_mass;
}

template <unsigned int dim>
void ParticleCluster<dim>::garbage_collect()
{
	/* We will now check whether we could optimise our particles.
		If some of our sub-clusters have particles, but less than the maximum among all of them, we can take those particles and
			transfer them to us.
		We only check for this condition after a certain depth since doing it from the root would take too long and would
			also be useless.
	*/
	if (!contains_particles())
	{
		if (nest_depth >= (log2(total_particles) / log2(max_children_particles)) && get_children_particle_num_recursive() < max_children_particles)
		{
			__gather_particles_recursive();
			active = true;
			has_particles = true;

			assert(children_clusters.size() == 0);
			assert(children_particles.size() <= max_children_particles);
		}
	}
	
#pragma omp parallel for
	for (long i = 0; i < children_clusters.size(); ++i)
	{
		auto c = std::next(children_clusters.begin(), i);
		if (!c->second.is_active())
		{
			c->second.remove_all_subclusters_recursive();
		}
		else
		{
			c->second.garbage_collect();
		}
	}

		for (auto it = children_clusters.begin(); it != children_clusters.end();)
		{
			if (!it->second.is_active())
				it = children_clusters.erase(it);
			else
				++it;
		}
	
}