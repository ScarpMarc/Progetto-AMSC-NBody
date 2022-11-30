#include <iostream>
#include <array>
#include <vector>
#include <memory>
#include "Particle.h"
#include "Vector.h"

#define DIM 3

using namespace std;

int main()
{

	int total_particles(3);

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
		for (unsigned int i = 0; i < DIM; i++)
		{
			std::cout << particles[i] ->get_position()[i] << std::endl;
		}

		std::cout << "With velocity" << std::endl;
		for (unsigned int i = 0; i < DIM; i++)
		{
			std::cout << particles[i] ->get_speed()[i] << std::endl;
		}

		std::cout << "and acceleration" << std::endl;
		for (unsigned int i = 0; i < DIM; i++)
		{
			std::cout << particles[i] ->get_acc()[i] << std::endl;
		}
	}
	
	
}