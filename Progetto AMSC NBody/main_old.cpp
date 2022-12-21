#include <iostream>
#include <array>
#include <vector>
#include <memory>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include "Particle.h"
#include "Vector.h"
#include "ForceMatrix.h"
#include "Constants.h"
#include "simulation_functions.h"
#include "json_parser.hpp"
#include "Profiling.h"
#include "Constants.h"

#include <random>
#include <numeric>
#include <thread>

#include "OpenGLFunctions.h"

#include "Globals.h"

using namespace std;

void saveParticles(const std::vector<Particle<DIM>>&, const std::string&);

void loadParticles(std::vector<Particle<DIM>>&, const std::string&);


// genera test case e salvalo su file di output
int main(int argc, char ** argv)
{
	std::string filename = "particles.pt";

	JsonParser parser("");
    parser.parse();

	// vector of Particle objects
	std::vector<Particle<DIM>> particles;
	Vector<DIM> position, speed, acceleration;


	for (int i = 0; i < total_particles; i++)
	{
		for (int j = 0; j < total_particles; j++)
		{
			for (int k = 0; k < total_particles; k++)
			{
				// generate mass
				// 1/9 3/9 5/9 7/9 9/9
				double mass(static_cast<double>(1.0e15));
				// generate new position, velocity and acceleration
				position = Vector<DIM>({ 2000.0 * ((double)i * 2.0 - 9.0) / 9.0,2000.0 * ((double)j * 2.0 - 9.0) / 9.0 ,2000.0 * ((double)k * 2.0 - 9.0) / 9.0 });
				speed = Vector<DIM>({ 0,0,0 });
				acceleration = Vector<DIM>({ 0,0,0 });

				// generate particle
				particles.push_back(<Particle<DIM>>(i, position, speed, acceleration, mass));
			}
		}
	}

	saveParticles(particles, filename);

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

void loadParticles(std::vector<Particle<DIM>> & particles, const std::string & filename)
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
