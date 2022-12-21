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

#include <random>
#include <numeric>
#include <thread>


#include "Globals.h"

using namespace std;

void saveParticles(const std::vector<Particle<DIM>>&, const std::string&);

void loadParticles(std::vector<Particle<DIM>>&, const std::string&);

unsigned int total_particles = 1000;
unsigned long int max_ticks = 100;
double ticks_per_second = 10.0;
unsigned int save_status_interval = 10;
unsigned int screen_refresh_millis = 200;
unsigned int screenResX = 2048;
unsigned int screenResY = 2048;

std::string save_filename = "particles_output.pt";
std::string load_filename;


// genera test case e salvalo su file di output
int main(int argc, char ** argv)
{
	std::string filename = "particles.pt";

	JsonParser parser("");
    parser.parse();

	// vector of Particle objects
	std::vector<Particle<DIM>> particles;
	Vector<DIM> position, speed, acceleration;

	int cuberoot_of_total_particles = 10;


	for (int i = 0; i < cuberoot_of_total_particles; i++)
	{
		for (int j = 0; j < cuberoot_of_total_particles; j++)
		{
			for (int k = 0; k < cuberoot_of_total_particles; k++)
			{
				// generate mass
				// 1/9 3/9 5/9 7/9 9/9
				double mass(static_cast<double>(1.0e15));
				// generate new position, velocity and acceleration
				position = Vector<DIM>({ 2000.0 * ((double)i * 2.0 - 9.0) / 9.0,2000.0 * ((double)j * 2.0 - 9.0) / 9.0 ,2000.0 * ((double)k * 2.0 - 9.0) / 9.0 });
				speed = Vector<DIM>({ 0,0,0 });
				acceleration = Vector<DIM>({ 0,0,0 });

				unsigned int ID = static_cast<unsigned int>(i);

				Particle<DIM> temp(ID, position,  speed, acceleration, mass);

				// generate particle
				particles.push_back(temp);
			}
		}
	}

	saveParticles(particles, filename);

}

void saveParticles(const std::vector<Particle<DIM>>& particles, const std::string& filename)
{
	std::ofstream outfile(filename, std::ios::out | std::ios::binary);
	if (!outfile)
	{
		std::cout << "Cannot open file to write particles!" << std::endl;
		return;
	}

	outfile.write(reinterpret_cast<const char*>(&DIM), sizeof(DIM));
	outfile.write(reinterpret_cast<char*>(&total_particles), sizeof(total_particles));

	for (const Particle<DIM>& particle : particles)
	{
		particle.saveToFile(outfile);
	}
	outfile.close();
}
// Close OpenGL window and terminate GLFW


void loadParticles(std::vector<Particle<DIM>>& particles, const std::string& filename)
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

	infile.read(reinterpret_cast<char*>(&total_particles), sizeof(unsigned int));
	for (int i = 0; i < total_particles; i++)
	{
		Particle<DIM> particle(i);
		particle.loadFromFile(infile);
		particles.emplace_back(particle);
	}
	infile.close();
}