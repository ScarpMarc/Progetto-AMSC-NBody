#include <iostream>
#include <array>
#include <vector>
#include <memory>
#include "Particle.h"
#include "Vector.h"
#include "ForceMatrix.h"
#include "Constants.h"

#include <thread>

//#define DELTA_T 0.1f
//#define MAX_TIME 10

using namespace std;


int main()
{
	const unsigned int total_particles(3);
	std::string filename = "particles.pt";

	// vector of unique pointers to Particle objects
	std::vector<std::unique_ptr<Particle<DIM>>> particles;
	Vector<DIM> position, speed, acceleration;

	for (unsigned int i = 0; i < total_particles; i++)
	{
		// generate mass
		double mass(static_cast<double>((/*i+*/1)*3*10e3));
		// generate new position, velocity and acceleration
		position = Vector<DIM>({ 0.0 + static_cast<double>(i),0.0 + static_cast<double>(i),0.0 + static_cast<double>(i) });
		//position = Vector<DIM>({ (i==1)*1.0 + (i == 2) * 0.7, (i==2)*0.5, 0.0});
		speed = Vector<DIM>({ 0.0,0.0,0.0 });
		acceleration = Vector<DIM>({ 0.0,0.0,0.0 });
		

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
	}
*/

	
	// UPDATE CYCLE

	ForceMatrix<DIM> force_matrix(total_particles);
	force_matrix.updateForces(particles);
	Vector<DIM> temp;
	unsigned int time(0);

	// UPDATE SECTION
	for (time = 0; time < max_ticks; ++time)
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
			particles[i]->calcNewPosition(1);
		}

		if(time == 0 || time == max_ticks - 1){
		
			cout << "----------------------------------------------------------------\nCycle " 
				<< time << " (Time " << time/ticks_per_second <<" seconds)\n" << endl;
			for (unsigned int i = 0; i < total_particles; i++)
			{
				std::cout << "Particle #:" << particles[i]->get_particle_id() << std::endl;
				std::cout << "In position" << std::endl;
				for (unsigned int j = 0; j < DIM; j++)
				{
					std::cout << particles[i]->get_position()[j] << std::endl;
				}

				std::cout << "With velocity" << std::endl;
				for (unsigned int j = 0; j < DIM; j++)
				{
					std::cout << particles[i]->get_speed()[j] << std::endl;
				}

				std::cout << "and acceleration" << std::endl;
				for (unsigned int j = 0; j < DIM; j++)
				{
					std::cout << particles[i]->get_acc()[j] << std::endl;
				}
			}
		}
	}	
}

void saveParticles(const std::vector<Particle<DIM>> & particles, const std::string & filename)
{
	std::ofstream outfile(filename, std::ios::out | std::ios::binary);
	if (!outfile)
	{
		std::cout << "Cannot open file to write particles!" << std::endl;
		return;
	}
	outfile << DIM;
	outfile << particles.size();
	for (Particle particle : particles)
	{
		particle.saveToFile(outfile);
	}
	outfile.close();
}

void laodParticles(std::vector<Particle<DIM>> & particles, const std::string & filename)
{
	std::ifstream infile(filename, std::ios::in | std::ios::binary);
	if (!infile)
	{
		std::cout << "Cannot open file to read particles!" << std::endl;
		return;
	}

	unsigned int dim;
	infile.read(reinterpret_cast<char*>(&dim), sizeof(unsigned int));
	if (dim != DIM)
	{
		std::cout << "File has a different particle DIM.\nCannot load file information!" << std::endl;
		return;
	}
	
	infile.read(reinterpret_cast<char*>(&dim), sizeof(unsigned int));
	for (int i = 0; i < dim; i++)
	{
		Particle<DIM> particle(i);
		particle.loadFromFile(infile);
		particles.emplace_back(particle);
	}
	infile.close();
}
