#pragma once
// #define GLEW_STATIC

#include <iostream>
#include <vector>

// #include <memory>

#include "Particle.h"
#include "ParticleCluster.h"
#include "shader.h"
#include "Constants.h"
#include "Controls.h"

#include <iostream>
#ifdef _WIN32
//#include <GL/gl.h>
#else

#endif

#include "../libs/include/GL/glew.h"
#include "../glfw/include/GLFW/glfw3.h"
#include "../glm/glm/glm.hpp"
#include "../glm/glm/gtc/matrix_transform.hpp"

#include <chrono>
#include <thread>
#include <array>
#include <random>
#include <cstring>

#include <iostream>

int gl_init(GLFWwindow** window);

extern "C"
{
	GLuint loadDDS(const char* imagepath);
}

/*template<unsigned int dim>
void create_cube_OLD(const Vector<dim>& m, const Vector<dim>& M, GLfloat* out_buffer_data)
{
	GLfloat Mx = static_cast<GLfloat>(M[0]) * 0.5 / screenResX;
	GLfloat My = static_cast<GLfloat>(M[1]) * 0.5 / screenResY;
	GLfloat Mz = static_cast<GLfloat>(M[2]) * 0.5 / screenResX;
	GLfloat mx = static_cast<GLfloat>(m[0]) * 0.5 / screenResX;
	GLfloat my = static_cast<GLfloat>(m[1]) * 0.5 / screenResY;
	GLfloat mz = static_cast<GLfloat>(m[2]) * 0.5 / screenResX;


	GLfloat temp_arr[] = {
	mx,my,mz,Mx,my,mz,mx,My,mz, // triangle 0 
	Mx,my,mz,mx,My,mz,Mx,My,mz, // 1
	Mx,my,mz,Mx,My,mz,Mx,my,Mz, // 2 
	Mx,My,mz,Mx,my,Mz,Mx,My,Mz, // 3 
	Mx,my,Mz,mx,my,Mz,mx,My,Mz, // 4
	Mx,my,Mz,Mx,My,Mz,mx,My,Mz, // 5
	mx,my,mz,mx,my,Mz,mx,My,mz, // 6
	mx,my,Mz,mx,My,mz,mx,My,Mz, // 7
	mx,my,mz,Mx,my,mz,Mx,my,Mz, // 8
	mx,my,mz,mx,my,Mz,Mx,my,Mz, // 9
	mx,My,mz,mx,My,Mz,Mx,My,Mz, // 10
	mx,My,mz,Mx,My,mz,Mx,My,Mz // 11
	};

	memcpy(out_buffer_data, temp_arr, sizeof(temp_arr));
}

template<unsigned int dim>
void create_cube(GLfloat* out_buffer_data)
{
	GLfloat Mx = 0.5f;
	GLfloat My = 0.5f;
	GLfloat Mz = 0.5f;
	GLfloat mx = -0.5f;
	GLfloat my = -0.5f;
	GLfloat mz = -0.5f;


	GLfloat temp_arr[] = {
	mx,my,mz,Mx,my,mz,mx,My,mz, // triangle 0 
	Mx,my,mz,mx,My,mz,Mx,My,mz, // 1
	Mx,my,mz,Mx,My,mz,Mx,my,Mz, // 2 
	Mx,My,mz,Mx,my,Mz,Mx,My,Mz, // 3 
	Mx,my,Mz,mx,my,Mz,mx,My,Mz, // 4
	Mx,my,Mz,Mx,My,Mz,mx,My,Mz, // 5
	mx,my,mz,mx,my,Mz,mx,My,mz, // 6
	mx,my,Mz,mx,My,mz,mx,My,Mz, // 7
	mx,my,mz,Mx,my,mz,Mx,my,Mz, // 8
	mx,my,mz,mx,my,Mz,Mx,my,Mz, // 9
	mx,My,mz,mx,My,Mz,Mx,My,Mz, // 10
	mx,My,mz,Mx,My,mz,Mx,My,Mz // 11
	};

	memcpy(out_buffer_data, temp_arr, sizeof(temp_arr));
}*/

/*template <unsigned int dim>
void populate_cluster_vertex_array(ParticleCluster<dim>* cluster, GLubyte* output, const unsigned int& index)
{
	for (unsigned int i = 0; i < dim; ++i)
	{
		output[index * 4 + i] = cluster ->get_global_max_boundary
	}
}*/

template <unsigned int dim>
void drawParticles(GLFWwindow** window, std::vector<Particle<dim>>* particles/*, ParticleCluster<dim>* main_cluster*/)
{
	std::random_device rd; // obtain a random number from hardware
	std::mt19937 gen(rd()); // seed the generator
	std::uniform_int_distribution<> distr(0, 6); // define the range

	std::array<std::array<unsigned char, 3>, 7> stars = {
		155, 176, 255,
		170, 191, 255,
		202, 215, 255,
		248, 247, 255,
		255, 244, 234,
		255, 210, 161,
		255, 204, 111 };

	GLuint VertexArrayID;
	glGenVertexArrays(1, &VertexArrayID);
	glBindVertexArray(VertexArrayID);

	// Create and compile our GLSL program from the shaders
	GLuint programID_main = LoadShaders("vertexShader.vertexshader", "fragmentShader.fragmentshader");
	//GLuint programID_simple = LoadShaders("simple_vertexshader.vertexshader", "simple_fragmentshader.fragmentshader");

	// Vertex shader
	GLuint CameraRight_worldspace_ID = glGetUniformLocation(programID_main, "CameraRight_worldspace");
	GLuint CameraUp_worldspace_ID = glGetUniformLocation(programID_main, "CameraUp_worldspace");
	GLuint ViewProjMatrixID = glGetUniformLocation(programID_main, "VP");

	// fragment shader
	GLuint TextureID = glGetUniformLocation(programID_main, "myTextureSampler");

	GLuint Texture = loadDDS("particle.DDS");

	static GLfloat* g_particule_position_size_data = new GLfloat[particles->size() * 4];
	static GLubyte* g_particule_color_data = new GLubyte[particles->size() * 4];
	//static GLubyte* g_cluster_color_data = new GLubyte[clusters->size() * 4];

	// The VBO containing the 4 vertices of the particles.
	// Thanks to instancing, they will be shared by all particles.
	static const GLfloat g_vertex_buffer_data[] = {
		-0.5f, -0.5f, 
		0.0f, 0.5f,
		-0.5f, 0.0f,
		-0.5f, 0.5f, 
		0.0f, 0.5f, 
		0.5f, 0.0f,
	};

	GLuint billboard_vertex_buffer;
	glGenBuffers(1, &billboard_vertex_buffer);
	glBindBuffer(GL_ARRAY_BUFFER, billboard_vertex_buffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(g_vertex_buffer_data), g_vertex_buffer_data, GL_STATIC_DRAW);

	// The VBO containing the positions and sizes of the particles
	GLuint particles_position_buffer;
	glGenBuffers(1, &particles_position_buffer);
	glBindBuffer(GL_ARRAY_BUFFER, particles_position_buffer);
	// Initialize with empty (NULL) buffer : it will be updated later, each frame.
	glBufferData(GL_ARRAY_BUFFER, particles->size() * 4 * sizeof(GLfloat), NULL, GL_STREAM_DRAW);

	// The VBO containing the colors of the particles
	GLuint particles_color_buffer;
	glGenBuffers(1, &particles_color_buffer);
	glBindBuffer(GL_ARRAY_BUFFER, particles_color_buffer);
	// Initialize with empty (NULL) buffer : it will be updated later, each frame.
	glBufferData(GL_ARRAY_BUFFER, particles->size() * 4 * 4 * sizeof(GLubyte), NULL, GL_STREAM_DRAW);

	// Generate second VAO
	/*GLuint VertexArrayID_boundaries;
	glGenVertexArrays(1, &VertexArrayID_boundaries);
	glBindVertexArray(VertexArrayID_boundaries);

	// Field boundaries
	GLfloat* boundary_vertices = (GLfloat*)malloc(12 * 3 * 3 * sizeof(GLfloat)); // 12 triangles, 3 vertices per triangle, 3 coords per vertex
	create_cube(Particle<dim>::get_global_max_boundary(), Particle<dim>::get_global_min_boundary(), boundary_vertices);

	// The VBO containing the geometry of a cluster (cube). This will not change, similarly to the billboard VBO for particles
	GLuint billboard_vertex_buffer_boundaries;
	glGenBuffers(1, &billboard_vertex_buffer_boundaries);
	glBindBuffer(GL_ARRAY_BUFFER, billboard_vertex_buffer_boundaries);
	// Initialize with empty (NULL) buffer : it will be updated later, each frame.
	glBufferData(GL_ARRAY_BUFFER, sizeof(boundary_vertices), boundary_vertices, GL_STATIC_DRAW);

	// The VBO containing the boundaries of the field. Will change based on which cluster we are drawing.
	GLuint global_boundaries_buffer;
	glGenBuffers(1, &global_boundaries_buffer);
	glBindBuffer(GL_ARRAY_BUFFER, global_boundaries_buffer);
	// Initialize with empty (NULL) buffer : it will be updated later, each frame.
	//glBufferData(GL_ARRAY_BUFFER, clusters->size() * 4 * 4 * sizeof(GLubyte), NULL, GL_STREAM_DRAW); 
	*/
	do
	{
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clear the screen

		// ------- Drawing particles
		glBindVertexArray(VertexArrayID);

		computeMatricesFromInputs(window);
		glm::mat4 ProjectionMatrix = getProjectionMatrix();
		glm::mat4 ViewMatrix = getViewMatrix();
		glm::vec3 CameraPosition(glm::inverse(ViewMatrix)[3]);
		glm::mat4 ViewProjectionMatrix = ProjectionMatrix * ViewMatrix;
		int ParticlesCount = 0;

		for (int i = 0; i < particles->size(); i++)
		{
			Particle<dim>& p = particles->at(i); // shortcut

			// p.cameradistance = glm::length2(p.pos - CameraPosition);

			// Fill the GPU buffer
			// TODO account for different dims
			g_particule_position_size_data[4 * i + 0] = (GLfloat)(p.get_position()[0]) / screenResX;
			g_particule_position_size_data[4 * i + 1] = (GLfloat)(p.get_position()[1]) / screenResY;
			g_particule_position_size_data[4 * i + 2] = (GLfloat)(p.get_position()[2]) / screenResX;
			g_particule_position_size_data[4 * i + 3] = .01; // p.getMass();

			//int thisColIdx = distr(gen);
			//g_particule_color_data[4 * i + 0] = stars[thisColIdx][0];	 // p.r;
			//g_particule_color_data[4 * i + 1] = stars[thisColIdx][1];	 // p.g;
			//g_particule_color_data[4 * i + 2] = stars[thisColIdx][2];	 // p.b;
			//g_particule_color_data[4 * i + 3] = 192; // p.a;
			g_particule_color_data[4 * i + 0] = 0;	 // p.r;
			g_particule_color_data[4 * i + 1] = 0;	 // p.g;
			g_particule_color_data[4 * i + 2] = 0;	 // p.b;
			g_particule_color_data[4 * i + 3] = 192; // p.a;

			++ParticlesCount;
		}

		glBindBuffer(GL_ARRAY_BUFFER, particles_position_buffer);
		glBufferData(GL_ARRAY_BUFFER, particles->size() * 4 * sizeof(GLfloat), NULL, GL_STREAM_DRAW); // Buffer orphaning, a common way to improve streaming perf.
		glBufferSubData(GL_ARRAY_BUFFER, 0, ParticlesCount * sizeof(GLfloat) * 4, g_particule_position_size_data);

		glBindBuffer(GL_ARRAY_BUFFER, particles_color_buffer);
		glBufferData(GL_ARRAY_BUFFER, particles->size() * 4 * sizeof(GLubyte), NULL, GL_STREAM_DRAW); // Buffer orphaning, a common way to improve streaming perf.
		glBufferSubData(GL_ARRAY_BUFFER, 0, ParticlesCount * sizeof(GLubyte) * 4, g_particule_color_data);

		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		// Use our shader
		glUseProgram(programID_main);

		// Bind our texture in Texture Unit 0
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, Texture);
		// Set our "myTextureSampler" sampler to use Texture Unit 0
		glUniform1i(TextureID, 0);

		// Same as the billboards tutorial
		glUniform3f(CameraRight_worldspace_ID, ViewMatrix[0][0], ViewMatrix[1][0], ViewMatrix[2][0]);
		glUniform3f(CameraUp_worldspace_ID, ViewMatrix[0][1], ViewMatrix[1][1], ViewMatrix[2][1]);

		glUniformMatrix4fv(ViewProjMatrixID, 1, GL_FALSE, &ViewProjectionMatrix[0][0]);

		// 1rst attribute buffer : vertices
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, billboard_vertex_buffer);
		glVertexAttribPointer(
			0,		  // attribute. No particular reason for 0, but must match the layout in the shader.
			3,		  // size
			GL_FLOAT, // type
			GL_FALSE, // normalized?
			0,		  // stride
			(void*)0 // array buffer offset
		);

		// 2nd attribute buffer : positions of particles' centers
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, particles_position_buffer);
		glVertexAttribPointer(
			1,		  // attribute. No particular reason for 1, but must match the layout in the shader.
			dim + 1,  // size : x + y + z + size => 4
			GL_FLOAT, // type
			GL_FALSE, // normalized?
			0,		  // stride
			(void*)0 // array buffer offset
		);

		// 3rd attribute buffer : particles' colors
		glEnableVertexAttribArray(2);
		glBindBuffer(GL_ARRAY_BUFFER, particles_color_buffer);
		glVertexAttribPointer(
			2,				  // attribute. No particular reason for 1, but must match the layout in the shader.
			4,				  // size : r + g + b + a => 4
			GL_UNSIGNED_BYTE, // type
			GL_TRUE,		  // normalized?    *** YES, this means that the unsigned char[4] will be accessible with a vec4 (floats) in the shader ***
			0,				  // stride
			(void*)0		  // array buffer offset
		);

		// These functions are specific to glDrawArrays*Instanced*.
		// The first parameter is the attribute buffer we're talking about.
		// The second parameter is the "rate at which generic vertex attributes advance when rendering multiple instances"
		// http://www.opengl.org/sdk/docs/man/xhtml/glVertexAttribDivisor.xml
		glVertexAttribDivisor(0, 0); // particles vertices : always reuse the same 4 vertices -> 0
		glVertexAttribDivisor(1, 1); // positions : one per quad (its center)                 -> 1
		glVertexAttribDivisor(2, 1); // color : one per quad                                  -> 1

		// Draw the particules !
		// This draws many times a small triangle_strip (which looks like a quad).
		// This is equivalent to :
		// for(i in ParticlesCount) : glDrawArrays(GL_TRIANGLE_STRIP, 0, 4),
		// but faster.

		glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, particles->size());

		glDisableVertexAttribArray(0);
		glDisableVertexAttribArray(1);
		glDisableVertexAttribArray(2);

		/*cluster_amt = main_cluster->get_subclusters_num_recursive();

		glBindVertexArray(VertexArrayID_boundaries);
		for (int i = 0; i < clusters->size(); i++)
		{

		}

		glBindBuffer(GL_ARRAY_BUFFER, global_boundaries_buffer);
		glBufferData(GL_ARRAY_BUFFER, cluster_amt * 4 * sizeof(GLfloat), NULL, GL_STREAM_DRAW); // Buffer orphaning, a common way to improve streaming perf.
		glUseProgram(programID_simple);

		// Draw boundaries
		// 4th attribute buffer : cube vertices
		glEnableVertexAttribArray(1);

		glBindBuffer(GL_ARRAY_BUFFER, global_boundaries_buffer);
		glVertexAttribPointer(
			1,		  // attribute. No particular reason for 0, but must match the layout in the shader.
			3,		  // size
			GL_FLOAT, // type
			GL_FALSE, // normalized?
			0,		  // stride
			(void*)0 // array buffer offset
		);

		//glVertexAttribDivisor(3, 0);

		glDrawElements(GL_TRIANGLE_STRIP, 0, 6);
		glDisableVertexAttribArray(1);
		//glDisableVertexAttribArray(1);
		*/

		// Swap buffers
		glfwSwapBuffers(*window);
		glfwPollEvents();

	} // Check if the ESC key was pressed or the window was closed
	while (glfwGetKey(*window, GLFW_KEY_ESCAPE) != GLFW_PRESS &&
		glfwWindowShouldClose(*window) == 0);

	// Cleanup VBO and shader
	glDeleteBuffers(1, &particles_color_buffer);
	glDeleteBuffers(1, &particles_position_buffer);
	glDeleteBuffers(1, &billboard_vertex_buffer);
	//glDeleteBuffers(1, &global_boundaries_buffer);
	glDeleteProgram(programID_main);
	glDeleteTextures(1, &Texture);
	glDeleteVertexArrays(1, &VertexArrayID);
	//glDeleteVertexArrays(1, &VertexArrayID_boundaries);

	//free(boundary_vertices);
}