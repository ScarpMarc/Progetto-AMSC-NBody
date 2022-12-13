#include <iostream>
#include <array>
#include <vector>
#include <memory>
#include <chrono>
#include "Particle.h"
#include "Vector.h"
#include "ForceMatrix.h"
#include "Constants.h"

#include <random>

#include <thread>

#include "OpenGLFunctions.h"

//#define DELTA_T 0.1f
//#define MAX_TIME 10

using namespace std;

std::vector<std::unique_ptr<Particle<DIM>>> particles;

const unsigned int testpnum = 10;
void mainLoop()
{
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
				double mass(static_cast<double>(1.0e16));
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


	for (unsigned int i = 0; i < total_particles; i++)
	{
		//(*(particles[i])).print();
	}



	/*
		std::cout << "Particle #:" << particles[i]->get_particle_id() << std::endl;
		std::cout << "In position" << std::endl;
		for (unsigned int j = 0; j < DIM; j++)
		{
			std::cout << particles[i] ->get_position()[j] << std::endl;
		}

		std::cout << "With velocity" << std::endl;
		for (unsigned int j = 0; j < DIM; j++)
		{
			std::cout << particles[i] ->get_speed()[j] << std::endl;
		}

		std::cout << "and acceleration" << std::endl;
		for (unsigned int j = 0; j < DIM; j++)
		{
			std::cout << particles[i] ->get_acc()[j] << std::endl;
		}
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
		auto start = chrono::high_resolution_clock::now();

		cout << "Iteration " << std::setw(6) << time << " (simulation seconds: " << std::setw(4) << (double)time / ticks_per_second << ") --- ";
		//compute forces
		force_matrix.updateForces(particles);
#pragma omp parallel for
		for (unsigned int i = 0; i < total_particles; i++)
		{
			// updating forces
			temp = force_matrix.getTotalForceOnParticle(i);
			particles[i]->updateResultingForce(temp);
		}

#pragma omp parallel for
		for (unsigned int i = 0; i < total_particles; i++)
		{
			// updating positions
			particles[i]->calcNewPosition(1);
		}
		for (unsigned int i = 0; i < total_particles; i++)
		{

			//(*(particles[i])).print();
		}

		//std::this_thread::sleep_for(std::chrono::milliseconds(10));
		auto end = chrono::high_resolution_clock::now();
		auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
		cout << "Execution time: " << std::setw(15) << duration.count() << " us" << endl;

	}

	auto simend = chrono::high_resolution_clock::now();
	auto simduration = chrono::duration_cast<chrono::microseconds>(simend - simstart);

	cout << "SIMULATION ENDED. Time taken: " << simduration.count() / 1000000 << " s" << endl;
}

int main()
{
	GLFWwindow* window = nullptr;
	gl_init(&window);

	std::thread t0(mainLoop);

	//drawParticles(&window, &particles);
	//std::thread t1(&drawParticles<DIM>, &window, &particles);

	// Close OpenGL window and terminate GLFW
	glfwTerminate();

	// TODO
	// Stop the other thread gracefully...
	t0.join();
}