#include <iostream>
#include <array>
#include <vector>
#include <memory>
#include "Particle.h"
#include "Vector.h"
#include "ForceMatrix.h"

#define DIM 3
#define DELTA_T 0.1
#define MAX_TIME 10

using namespace std;


int main()
{

	const unsigned int total_particles(3);

	// vector of unique pointers to Particle objects
	std::vector<std::unique_ptr<Particle<DIM>>> particles;
	double mass;
	Vector<DIM> position, speed, acceleration;

	for (unsigned int i = 0; i < total_particles; i++)
	{

		// generate mass
		double mass(static_cast<double>(i*3));
		// generate new position, velocity and acceleration
		position = Vector<DIM>({0.0 + static_cast<double>(i),0.0 + static_cast<double>(i),0.0 + static_cast<double>(i)});
		speed = Vector<DIM>({0.0,0.0,0.0});
		acceleration = Vector<DIM>({0.0,0.0,0.0});

		// generate particle
		particles.push_back(std::make_unique<Particle<DIM>>(i, position, speed, acceleration, mass));
	}
	
	// print particles

	for (unsigned int i = 0; i < total_particles; i++)
	{
		std::cout << "Particle #:" << particles[i] ->get_particle_id() << std::endl;
		std::cout << "In position" << std::endl;
		for (unsigned int j = 0; j < DIM; j++)
		{
			std::cout << particles[i] ->get_position()[j] << std::endl;
		}

		std::cout << "With velocity" << std::endl;
		for (unsigned int j = 0; j < DIM; j++)
		{
			std::cout << particles[j] ->get_speed()[j] << std::endl;
		}

		std::cout << "and acceleration" << std::endl;
		for (unsigned int j = 0; j < DIM; j++)
		{
			std::cout << particles[j] ->get_acc()[j] << std::endl;
		}
	}



	// UPDATE CYCLE

	ForceMatrix<DIM> force_matrix = ForceMatrix<DIM>(total_particles);
	force_matrix.updateForces(particles);
	Vector<DIM> temp;
	unsigned int time(0);

	// UPDATE SECTION
	for (time = 0; time < MAX_TIME; time++)
	{
		//compute forces
		force_matrix.updateForces(particles);
		#pragma omp parallel for
		for (unsigned int i = 0; i < total_particles; i++)
		{
			// updating forces
			temp = force_matrix.getTotalForceOnParticle(i);
			particles[i] -> updateResultingForce(temp);
		}

		#pragma omp parallel for
		for (unsigned int i = 0; i < total_particles; i++)
		{
			// updating positions
			particles[i]->calcNewPosition(DELTA_T);
		}
		

		time += DELTA_T;
	}
	

	
	
}