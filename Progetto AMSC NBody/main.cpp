#include <iostream>
#include <array>
#include <vector>
#include <memory>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <algorithm>
#include "Particle.h"
#include "ParticleCluster.h"
#include "Vector.h"
#include "ForceMatrix.h"
#include "Constants.h"
#include "simulation_functions.h"
#include "json_parser.hpp"
#include "Profiling.h"

#include <random>
#include <numeric>
#include <thread>

#ifdef USE_GRAPHICS
#include "OpenGLFunctions.h"
#endif

#include "Globals.h"

// Must do initialisation in a source file
template<unsigned int dim>
unsigned int Particle<dim>::maxID = 0;

// Must do initialisation in a source file
template<unsigned int dim>
unsigned int ParticleCluster<dim>::maxID = 0;
template<unsigned int dim>
unsigned int ParticleCluster<dim>::max_children = 8;

template<unsigned int dim>
Vector<dim> Particle<dim>::max_boundary; 
template<unsigned int dim>
Vector<dim> Particle<dim>::min_boundary;

using namespace std;

long long int forceComp_mean_durations_per_tick = 0, posComp_mean_durations_per_tick = 0, matrixComp_mean_duration = 0;
void saveParticles(const std::vector<Particle<DIM>> &, const std::string &);

extern long long int total_sim_duration = 0;
time_t programme_start;
std::string profiling_folder = "";
std::string profiling_file_name = "Profiler_.txt";

// vector of Particle objects
std::vector<Particle<DIM>> particles;

void loadParticles(std::vector<Particle<DIM>> &, const std::string &);

unsigned int total_particles = 1000;
unsigned long int max_ticks = 100;
double ticks_per_second = 10.0;
unsigned int save_status_interval = 10;
unsigned int screen_refresh_millis = 200;
unsigned int screenResX = 2048;
unsigned int screenResY = 2048;
bool use_graphics = true;
double domain_start_x = -5000.0, domain_end_x = 5000.0;
double domain_start_y = -5000.0, domain_end_y = 5000.0;
double domain_start_z = -5000.0, domain_end_z = 5000.0;

std::string save_filename = "particles_output.pt";
std::string load_filename;

bool load_particles_from_file = false;

#ifdef _WIN32
extern "C"
{
	__declspec(dllexport) unsigned long NvOptimusEnablement = 0x00000001;
	__declspec(dllexport) int AmdPowerXpressRequestHighPerformance = 1;
}
#endif

int mainLoop()
{
	total_sim_duration = 0;

	if (load_particles_from_file)
	{
		std::cout << "Attempting to load particle file..." << std::endl;
		loadParticles(particles, load_filename);
		std::cout << "Particle file loaded." << std::endl;
	}
	else
	{
		std::cout << "Attempting to load JSON file..." << std::endl;
		JsonParser parser("");
		parser.parse();

		std::srand(std::time(NULL));

		programme_start = time(0);

		std::cout << "JSON file loaded." << std::endl;

		for (int i = 0; i < total_particles; i++)
		{
			Vector<DIM> position, speed, acceleration;
			// generate mass
			double mass(static_cast<double>((std::rand() % 10000) * 1.0e7 + 1.0e7));
			// generate new position, velocity and acceleration
			position = Vector<DIM>({static_cast<double>(std::rand() % 4000) - 2000.0, static_cast<double>(std::rand() % 4000) - 2000.0, static_cast<double>(std::rand() % 4000) - 2000.0});
			// position = Vector<DIM>({ (i==1)*1.0 + (i == 2) * 0.7, (i==2)*0.5, 0.0});
			speed = Vector<DIM>({0.0, 0.0, 0.0});
			acceleration = Vector<DIM>({0.0, 0.0, 0.0});

			// generate particle
			particles.emplace_back(position, speed, acceleration, mass);
		}

		std::cout << "Particles generated." << std::endl;
	}

	Vector<DIM> temp;
	unsigned int time(0);

	auto simstart = chrono::high_resolution_clock::now();

	// UPDATE SECTION
	for (time = 0; time < max_ticks; ++time)
	{
		if (!(time % (max_ticks / 100)))
		{
			cout << "Iteration " << std::setw(6) << time << " (simulation seconds: " << std::setw(4) << (double)(time + 1) / ticks_per_second << ")";
		}

		std::chrono::microseconds matrixComp_duration_this_tick;
		auto matrixComp_start = chrono::high_resolution_clock::now();

		do_simulation_step(particles, 1);

		// compute forces
		auto matrixComp_end = chrono::high_resolution_clock::now();
		matrixComp_duration_this_tick = chrono::duration_cast<chrono::microseconds>(matrixComp_end - matrixComp_start);
		matrixComp_mean_duration += matrixComp_duration_this_tick.count();

		if (!(time % (max_ticks / 100)))
		{
			cout << " --- Execution time: " << std::setw(15) << matrixComp_duration_this_tick.count() << " us";
			cout << endl;
		}

		if (!(time % save_status_interval))
		{
			saveParticles(particles, save_filename);
		}
		/*
		if (time == 0 || time == max_ticks - 1)
		{
			cout << "----------------------------------------------------------------\nCycle "
				<< time << " (Time " << time / ticks_per_second << " seconds)\n" << endl;
			for (unsigned int i = 0; i < total_particles; i++)
			{
				std::cout << "Particle #:" << particles[i].get_particle_id() << std::endl;
				std::cout << "In position" << std::endl;
				for (unsigned int j = 0; j < DIM; j++)
				{
					std::cout << particles[i].get_position()[j] << std::endl;
				}
				std::cout << "With velocity" << std::endl;
				for (unsigned int j = 0; j < DIM; j++)
				{
					std::cout << particles[i].get_speed()[j] << std::endl;
				}

				std::cout << "and acceleration" << std::endl;
				for (unsigned int j = 0; j < DIM; j++)
				{
					std::cout << particles[i].get_acc()[j] << std::endl;
				}
			}
		}*/

	} // For time
	// Calculate mean time
	matrixComp_mean_duration /= max_ticks;

	// Calculate total execution time
	auto simend = chrono::high_resolution_clock::now();
	auto simduration = chrono::duration_cast<chrono::microseconds>(simend - simstart);

	cout << "SIMULATION ENDED. Time taken: " << (double)simduration.count() / 1000000.0 << " s" << endl;
	total_sim_duration = simduration.count();

	save_profiler_data_text_file(profiling_folder + profiling_file_name);

	return 0;
}

int main(int argc, char **argv)
{
#ifdef USE_GRAPHICS
	use_graphics = true;
#else
	use_graphics = false;
#endif
	/*if (argc != 3 && argc != 2)
	{
		cout << "Insufficient number of parameters!" << endl;
		return 1;
	}

	if (string(argv[1]) == "-json")
	{
		load_particles_from_file = false;
	}
	else if (string(argv[1]) == "-load")
	{
		load_filename = string(argv[2]);
		load_particles_from_file = true;
	}
	else
	{
		cout << "First parameter not recognised!" << endl;
		return 2;
	}*/
#ifdef USE_GRAPHICS
	GLFWwindow* window = nullptr;
	if (use_graphics)
	{
		gl_init(&window);
	}
#endif

	if (!use_graphics)
	{
		cout << "Computing without graphics" << endl;
	}

	std::thread t0(mainLoop);

#ifdef USE_GRAPHICS
	if (use_graphics)
	{
		drawParticles(&window, &particles);
		// std::thread t1(&drawParticles<DIM>, &window, &particles);

		glfwTerminate();
	}
#endif

	// TODO
	// Stop the other thread gracefully...
	t0.join();
}

void saveParticles(const std::vector<Particle<DIM>> &particles, const std::string &filename)
{
	std::ofstream outfile(filename, std::ios::out | std::ios::binary);
	if (!outfile)
	{
		std::cout << "Cannot open file to write particles!" << std::endl;
		return;
	}

	outfile.write(reinterpret_cast<const char *>(&DIM), sizeof(DIM));
	outfile.write(reinterpret_cast<char *>(&total_particles), sizeof(total_particles));

	for (const Particle<DIM> &particle : particles)
	{
		particle.saveToFile(outfile);
	}
	outfile.close();
}
// Close OpenGL window and terminate GLFW

void loadParticles(std::vector<Particle<DIM>> &particles, const std::string &filename)
{
	std::ifstream infile(filename, std::ios::in | std::ios::binary);
	if (!infile)
	{
		std::cout << "Cannot open file to read particles!" << std::endl;
		return;
	}

	unsigned int dim;
	infile.read(reinterpret_cast<char *>(&dim), sizeof(unsigned int));
	if (dim != DIM)
	{
		std::cout << "File has a different particle DIM.\nCannot load file information!" << std::endl;
		return;
	}

	infile.read(reinterpret_cast<char *>(&total_particles), sizeof(unsigned int));
	for (int i = 0; i < total_particles; i++)
	{
		Particle<DIM> particle;
		particle.loadFromFile(infile);
		particles.emplace_back(particle);
	}
	infile.close();
}
