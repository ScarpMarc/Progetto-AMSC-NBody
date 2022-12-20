#include <iostream>
#include <array>
#include <vector>
#include <memory>
#include <chrono>
#include "Particle.h"
#include "Vector.h"
#include "ForceMatrix.h"
#include "Constants.h"
#include "Profiling.h"

#include <random>
#include <numeric>
#include <thread>

#include "OpenGLFunctions.h"

#include "Globals.h"

using namespace std;

#ifdef ADVANCED_PROFILING
extern long long int forceComp_mean_durations_per_tick = 0, posComp_mean_durations_per_tick = 0, matrixComp_mean_duration = 0;

#endif
extern long long int total_sim_duration = 0;
time_t programme_start;
std::string profiling_folder = "";
std::string profiling_file_name = "Profiler_.txt";

std::vector<std::unique_ptr<Particle<DIM>>> particles;

const unsigned int testpnum = 10;
void mainLoop()
{
	total_sim_duration = 0;

	programme_start = time(0);
	// vector of unique pointers to Particle objects
	Vector<DIM> position, speed, acceleration;

	std::random_device rd; // obtain a random number from hardware
	std::mt19937 gen(rd()); // seed the generator
	std::uniform_int_distribution<> distr(-2000, 2000); // define the range

	for (int i = 0; i < testpnum; i++)
	{
		for (int j = 0; j < testpnum; j++)
		{
			for (int k = 0; k < testpnum; k++)
			{
				// generate mass
				// 1/9 3/9 5/9 7/9 9/9
				double mass(static_cast<double>(1.0e15));
				// generate new position, velocity and acceleration
				position = Vector<DIM>({ 2000.0 * ((double)i * 2.0 - 9.0) / 9.0,2000.0 * ((double)j * 2.0 - 9.0) / 9.0 ,2000.0 * ((double)k * 2.0 - 9.0) / 9.0 });
				speed = Vector<DIM>({ 0,0,0 });
				acceleration = Vector<DIM>({ 0,0,0 });

				// generate particle
				particles.push_back(std::make_unique<Particle<DIM>>(i, position, speed, acceleration, mass));
			}
		}
	}

	// print particles


	/*for (unsigned int i = 0; i < total_particles; i++)
	{
		//(*(particles[i])).print();
	}*/

	// UPDATE CYCLE

	ForceMatrix<DIM> force_matrix(total_particles);
	force_matrix.updateForces(particles);
	Vector<DIM> temp;
	unsigned int time(0);

	auto simstart = chrono::high_resolution_clock::now();

	// UPDATE SECTION
	for (time = 0; time < max_ticks; ++time)
	{
#ifdef ADVANCED_PROFILING
		std::array<long long int, total_particles> forceComp_durations_this_tick, posComp_durations_this_tick;
		std::chrono::microseconds matrixComp_duration_this_tick;
#endif
		cout << "Iteration " << std::setw(6) << time << " (simulation seconds: " << std::setw(4) << (double)(time + 1) / ticks_per_second << ")";

		auto start = chrono::high_resolution_clock::now();

#ifdef ADVANCED_PROFILING
		auto matrixComp_start = chrono::high_resolution_clock::now();
#endif
		//compute forces
		force_matrix.updateForces(particles);
#ifdef ADVANCED_PROFILING
		auto matrixComp_end = chrono::high_resolution_clock::now();
		matrixComp_duration_this_tick = chrono::duration_cast<chrono::microseconds>(matrixComp_end - matrixComp_start);
		matrixComp_mean_duration += matrixComp_duration_this_tick.count();
#endif

#pragma omp parallel for private(temp)
		for (int i = 0; i < total_particles; i++)
		{
#ifdef ADVANCED_PROFILING
			auto forceComp_start = chrono::high_resolution_clock::now();
#endif
			// updating forces
			temp = force_matrix.getTotalForceOnParticle(i);
			particles[i]->updateResultingForce(temp);

#ifdef ADVANCED_PROFILING
			auto forceComp_end = chrono::high_resolution_clock::now();
			forceComp_durations_this_tick[i] = chrono::duration_cast<chrono::microseconds>(forceComp_end - forceComp_start).count();
#endif
		}

#pragma omp parallel for
		for (int i = 0; i < total_particles; i++)
		{
#ifdef ADVANCED_PROFILING
			auto posComp_start = chrono::high_resolution_clock::now();
#endif
			// updating positions
			particles[i]->calcNewPosition(1);

#ifdef ADVANCED_PROFILING
			auto posComp_end = chrono::high_resolution_clock::now();
			posComp_durations_this_tick[i] = chrono::duration_cast<chrono::microseconds>(posComp_end - posComp_start).count();
#endif
		}

		/*for (unsigned int i = 0; i < total_particles; i++)
		{
			//(*(particles[i])).print();
		}*/

		forceComp_mean_durations_per_tick += accumulate(forceComp_durations_this_tick.begin(), forceComp_durations_this_tick.end(), 0LL) / total_particles;

		posComp_mean_durations_per_tick += accumulate(posComp_durations_this_tick.begin(), posComp_durations_this_tick.end(), 0LL) / total_particles;

		//std::this_thread::sleep_for(std::chrono::milliseconds(10));
		auto end = chrono::high_resolution_clock::now();
		auto duration = chrono::duration_cast<chrono::microseconds>(end - start);

		cout << " --- Execution time: " << std::setw(15) << duration.count() << " us";

		cout << endl;
	}

	auto simend = chrono::high_resolution_clock::now();
	auto simduration = chrono::duration_cast<chrono::microseconds>(simend - simstart);

	forceComp_mean_durations_per_tick /= max_ticks;
	posComp_mean_durations_per_tick /= max_ticks;
	matrixComp_mean_duration /= max_ticks;

	cout << "SIMULATION ENDED. Time taken: " << simduration.count() / 1000000 << " s" << endl;
	total_sim_duration = simduration.count();

	save_profiler_data_text_file(profiling_folder + profiling_file_name);
}

int main()
{
	GLFWwindow* window = nullptr;
	gl_init(&window);

	std::thread t0(mainLoop);

	drawParticles(&window, &particles);
	//std::thread t1(&drawParticles<DIM>, &window, &particles);

	// Close OpenGL window and terminate GLFW
	glfwTerminate();

	// TODO
	// Stop the other thread gracefully...
	t0.join();
}