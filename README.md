# Progetto-AMSC-NBody

# Description
This programme implements an N-Body problem solver. At the moment, the evaluation algorithm is the simple one with a complexity of O(n^2). The task for this project was to implement the FMM algorithm, but at present it is not finished; instead, I have managed to complete a smart clustering system for the particles. 

Particles are automatically divided in clusters that encompass the entire simulation field, which has no boundaries. Each cluster may either contain a fixed number of particles (currently 8, unfortunately not parametric as I did not have the time, but can be easily changed by modifying the appropriate variable in `ParticleCluster.h`) or be the root of a tree of sub-clusters (arranged in a cube, 2 per side). When particles are added to a cluster, it will check whether it can store them itself or "send them down" to one of its children, possibly spawning children on the fly. Evaluation of which children will be the recipient is done in O(1) by simply dividing the space in quadrants and indexing the children with respect to their relative position on the x, y (and z, depending on the problem dimension) axes.

At each frame, the main cluster resizes itself as a cube containing all the particles, and resizes all of its children recursively. Particles that "fall out" of clusters during this process, or "escape" while moving, are automatically redistributed by each cluster by first attempting to send them to their parent, which will evaluate whether to send it back down to another child or forward it to its own parent.

After all particles have been redistributed, all clusters calculate their center of mass, their total velocity (which is the vector sum of all velocities of sub-clusters or particles), etc. in a parallel fashion. This should enable the correct execution of the FMM algorithm in the future.

Every hundredth frame, a garbage-collector is called. It will traverse the cluster tree top-to-bottom and eliminate clusters that are not active (i.e. that have no particles within them, and also no particles within their sub-cluster tree); it will also check, past a certain "depth", whether the total count of particles within a cluster's tree of sub-cluser could be aggregated into itself. This ensures that not too many clusters are created, and also that there are no clusters that have "too few" particles among them.

# How to run
## Pulling
Always pull with `--recurse-submodules`. This programme needs `glew`, `glfw` and `glm` to be run with graphics.

## Required files
`settings.json` must be in the same directory as the executable. If compiled with graphics,
`vertexShader.vertexshader`, `fragmentShader.fragmentshader` and `particle.dds` must also be within the same directory.
These can be automatically copied to the right folder with `make deploy`

## JSON file options
- `DIM` (integer): problem dimension (2D or 3D) - UNTESTED
- `max_ticks` (integer): amount of ticks to run the simulation for
- `ticks_per_second` (floating-point): amount of ticks in a second of simulation. Values below 1 are possible (but always greater than 0); doing so will cause the simulation to "blow up"
- `screen_refresh_millis` (integer): screen refresh sleep time
- `screenResX` (integer): screen horizontal resolution
- `screenResY` (integer): screen vertical resolution
- `save_status_interval' (integer): delay in ticks between consecutive status saves
- `total_particles` (integer): amount of particles in the simulation

## Controls
Left and right arrows to strafe left and right; up and down arrows to zoom; ESC to close the graphics window. The rest of the programme should continue simulating, but there will be no way to re-open the window without restarting.

# How to build 
Skip sections 1 and 2 if you _do not_ want graphics. Alternatively, it is possible to inhibit graphics by specifying the CMake flag `-DNO_GRAPHICS=TRUE`.

**WARNING**: if you have installed the [mk-modules](https://github.com/elauksap/mk) and want to use the graphics, you **MUST** unload `gcc-glibc` as this interferes with the CMake routines that find OpenGL. In general, you should always refer to the latest version of CMake.
1. Go to `glew/auto` and `make`; this will configure the appropriate dependencies. You may need to access `Makefile` and `modify PYTHON ?= python` into `PYTHON ?= python3`
2. You may need to install the following libraries: `libxrandr-dev libxinerama-dev libxcursor-dev libxi-dev`
3. Create a `build` directory in the root folder and `cmake ..` inside of it. The output will state whether OpenGL has been found (if not, read the warning and refer to points 1 and 2).
4. `make`. The result should be in the `Progetto AMSC NBody` folder in `build`.
5. `make deploy` to automatically copy the required files in the right folder.
