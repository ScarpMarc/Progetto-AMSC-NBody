#pragma once
#include <iostream>
#include <vector>

#include "Particle.h"
#include "shader.h"

#include "../glew/Include/GL/glew.h"
#include "../glfw/include/GLFW/glfw3.h"
#include "../glm/glm/glm.hpp"

const unsigned int screenResX = 1024;
const unsigned int screenResY = 1024;

//const unsigned int dimX = 100;
//const unsigned int dimY = 100;

int gl_init(GLFWwindow* window);

template <unsigned int dim>
void drawParticles(GLFWwindow* window, const std::vector<Particle<dim>>& particles)
{
	GLuint VertexArrayID;
	glGenVertexArrays(1, &VertexArrayID);
	glBindVertexArray(VertexArrayID);
	// Create and compile our GLSL program from the shaders
	GLuint programID = LoadShaders("vertexShader.vertexshader", "fragmentShader.fragmentshader");

	//Position vector
	std::vector<GLfloat> vertexBufferData = std::vector<GLfloat>(particles.size() * dim, 0);
	// Stroke weight vector
	std::vector<GLfloat> weightVector = std::vector<GLfloat>(particles.size());

	// This will identify our vertex buffer
	GLuint vertexbuffer;
	// Generate 1 buffer, put the resulting identifier in vertexbuffer
	glGenBuffers(1, &vertexbuffer);
	// The following commands will talk about our 'vertexbuffer' buffer
	glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
	// Give our vertices to OpenGL.
	glBufferData(GL_ARRAY_BUFFER, vertexBufferData.size(), vertexBufferData.begin()._Ptr, GL_STATIC_DRAW);

	unsigned int iteration = 0;

	// TODO: proof of concept. Does not support adding new particles. Also: create external function for updating weights 
	//		and positions separately since weights are not expected to change often

	do
	{
		// Iterate over all particles and save each position and weight component
		{
			unsigned int k = 0, l = 0;
			for (unsigned int i = 0; i < particles.size(); ++i)
			{
				// Iterate over the dimensions
				for (unsigned int vectorDim = 0; vectorDim < dim; ++vectorDim)
				{
					vertexBufferData[k++] = (particles[i]).get_position()[vectorDim];

					weightVector[l++] = (particles[i]).get_mass();
				}
			}
		}

		// Clear the screen
		glClear(GL_COLOR_BUFFER_BIT);

		// Use our shader
		glUseProgram(programID);

		// 1st attribute buffer : vertices
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
		glVertexAttribPointer(
			0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
			dim,                // size
			GL_FLOAT,           // type
			GL_FALSE,           // normalized?
			0,                  // stride
			(void*)0            // array buffer offset
		);

		// Draw
		for (unsigned int counter = 0; counter < particles.size(); ++counter)
		{
			glPointSize((GLfloat)10 * weightVector[counter]);
			glDrawArrays(GL_POINTS, counter, 1);
		}

		glDisableVertexAttribArray(0);

		// Swap buffers
		glfwSwapBuffers(window);
		glfwPollEvents();

		iteration += 1;
		//Sleep(10);

	} // Check if the ESC key was pressed or the window was closed
	while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS &&
		glfwWindowShouldClose(window) == 0);

	// Cleanup VBO
	glDeleteBuffers(1, &vertexbuffer);
	glDeleteVertexArrays(1, &VertexArrayID);
	glDeleteProgram(programID);

	// Close OpenGL window and terminate GLFW
	glfwTerminate();

}