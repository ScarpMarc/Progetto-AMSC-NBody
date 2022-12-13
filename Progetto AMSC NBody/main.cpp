#include <iostream>
#include <array>
#include <vector>
#include <memory>
#include "Particle.h"
#include "Vector.h"
#include "ForceMatrix.h"
#include "Constants.h"

#include <random>

#include <thread>

//#define DELTA_T 0.1f
//#define MAX_TIME 10

using namespace std;

std::vector<std::unique_ptr<Particle<DIM>>> particles;


void mainLoop()
{
	// vector of unique pointers to Particle objects
	Vector<DIM> position, speed, acceleration;

	std::random_device rd; // obtain a random number from hardware
	std::mt19937 gen(rd()); // seed the generator
	std::uniform_int_distribution<> distr(1, 10); // define the range


	for (unsigned int i = 0; i < total_particles; i++)
	{
		// generate mass
		double mass(static_cast<double>((i + 1) * 10 * distr(gen)));
		// generate new position, velocity and acceleration
		position = Vector<DIM>({ 0.0 + static_cast<double>((double)distr(gen)),0.0 + static_cast<double>((double)distr(gen)),0.0 - static_cast<double>((double)distr(gen)) });
		speed = Vector<DIM>({ (double)distr(gen),(double)distr(gen),(double)distr(gen) });
		acceleration = Vector<DIM>({ (double)distr(gen),(double)distr(gen),(double)distr(gen) });

		// generate particle
		particles.push_back(std::make_unique<Particle<DIM>>(i, position, speed, acceleration, mass));
	}

	// print particles

	/*
	for (unsigned int i = 0; i < total_particles; i++)
	{
		(*(particles[i])).print();
	}
	*/


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

	// UPDATE SECTION
	for (time = 0; time < max_ticks; ++time)
	{
		cout << " --- Iteration " << time << " (simulation seconds: " << (double)time / ticks_per_second << ")" << endl;
		//compute forces
		force_matrix.updateForces(particles);
		//#pragma omp parallel for
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

		std::this_thread::sleep_for(std::chrono::milliseconds(10));

	}
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
}