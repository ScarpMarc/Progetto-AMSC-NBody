#pragma once
#include <iostream>
#include <vector>

//#include <memory>

#include "Particle.h"
#include "shader.h"

#include "../glew/Include/GL/glew.h"
#include "../glfw/include/GLFW/glfw3.h"
#include "../glm/glm/glm.hpp"
#include "../glm/glm/gtc/matrix_transform.hpp"

#include <Windows.h>

const unsigned int screenResX = 1024;
const unsigned int screenResY = 1024;

//const unsigned int dimX = 100;
//const unsigned int dimY = 100;

int gl_init(GLFWwindow** window);

template <unsigned int dim>
void drawParticles(GLFWwindow** window, const std::vector<Particle<dim>>& particles)
{
	GLuint VertexArrayID;
	glGenVertexArrays(1, &VertexArrayID);
	glBindVertexArray(VertexArrayID);
	// Create and compile our GLSL program from the shaders
	GLuint programID = LoadShaders("vertexShader.vertexshader", "fragmentShader.fragmentshader");

	//Position vector
	std::vector<GLfloat> vertexBufferData = std::vector<GLfloat>(particles.size() * dim);
	// Stroke weight vector
	std::vector<GLfloat> weightVector = std::vector<GLfloat>(particles.size());

	std::vector<GLfloat> colourVector = std::vector<GLfloat>(particles.size() * 3);
	{
		for (unsigned int i = 0; i < particles.size()*3; ++i)
		{
			colourVector[i] = ((double)rand() / (RAND_MAX));
		}
	}

	GLuint MatrixID = glGetUniformLocation(programID, "MVP");

	// Projection matrix : 45° Field of View, 4:3 ratio, display range : 0.1 unit <-> 10 units
	glm::mat4 Projection = glm::perspective(glm::radians(45.0f), 4.0f / 3.0f, .1f, 100.0f);
	// Camera matrix
	glm::mat4 View = glm::lookAt(
		glm::vec3(0, 0, 1), // Camera is at (4,3,-3), in World Space
		glm::vec3(0, 0, 0), // and looks at the origin
		glm::vec3(0, 1, 0)  // Head is up (set to 0,-1,0 to look upside-down)
	);
	// Model matrix : an identity matrix (model will be at the origin)
	glm::mat4 Model = glm::mat4(1.0f);
	// Our ModelViewProjection : multiplication of our 3 matrices
	glm::mat4 MVP = Projection * View * Model; // Remember, matrix multiplication is the other way around

	// This will identify our vertex buffer
	GLuint vertexbuffer;
	// Generate 1 buffer, put the resulting identifier in vertexbuffer
	glGenBuffers(1, &vertexbuffer);
	

	GLuint colourbuffer;
	glGenBuffers(1, &colourbuffer);
	glBindBuffer(GL_ARRAY_BUFFER, colourbuffer);
	glBufferData(GL_ARRAY_BUFFER, colourVector.size(), colourVector.begin()._Ptr, GL_STATIC_DRAW);

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
					vertexBufferData[k++] = (particles[i]).get_position()[vectorDim] / screenResX;
				}
				weightVector[l++] = (particles[i]).get_mass();
			}
		}

		// The following commands will talk about our 'vertexbuffer' buffer
		glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
		// Give our vertices to OpenGL.
		glBufferData(GL_ARRAY_BUFFER, vertexBufferData.size(), vertexBufferData.begin()._Ptr, GL_DYNAMIC_DRAW);

		// Clear the screen
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clear color and depth buffers
		//glMatrixMode(GL_MODELVIEW);     // To operate on model-view matrix

		// Send our transformation to the currently bound shader, 
		// in the "MVP" uniform
		glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);

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

		// 2nd attribute buffer : colors
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, colourbuffer);
		glVertexAttribPointer(
			1,                                // attribute. No particular reason for 1, but must match the layout in the shader.
			3,                                // size
			GL_FLOAT,                         // type
			GL_FALSE,                         // normalized?
			0,                                // stride
			(void*)0                          // array buffer offset
		);

		// Draw
		for (unsigned int counter = 0; counter < particles.size(); ++counter)
		{
			glPointSize((GLfloat) weightVector[counter]*2);
			glDrawArrays(GL_POINTS, counter, 1);
		}

		glDisableVertexAttribArray(0);
		glDisableVertexAttribArray(1);

		// Swap buffers
		glfwSwapBuffers(*window);
		glfwPollEvents();

		iteration += 1;
		Sleep(100);

	} // Check if the ESC key was pressed or the window was closed
	while (glfwGetKey(*window, GLFW_KEY_ESCAPE) != GLFW_PRESS &&
		glfwWindowShouldClose(*window) == 0);

	// Cleanup VBO
	glDeleteBuffers(1, &vertexbuffer);
	glDeleteVertexArrays(1, &VertexArrayID);
	glDeleteProgram(programID);

	glDeleteVertexArrays(1, &VertexArrayID);

	// Close OpenGL window and terminate GLFW
	glfwTerminate();

}