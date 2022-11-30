#include "OpenGLFunctions.h"
#include "shader.h"

int gl_init()
{
	TerrainField field = TerrainField(dimX, dimY);

	// GLFW initialisation
	if (!glfwInit())
	{
		return -1;
	}

	// Open window
	glfwWindowHint(GLFW_SAMPLES, 4); // 4x antialiasing
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4); // We want OpenGL 3.3
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy; should not be needed
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); // We don't want the old OpenGL 

	// Open a window and create its OpenGL context
	GLFWwindow* window; // (In the accompanying source code, this variable is global for simplicity)
	window = glfwCreateWindow(1024, 1024, "Landscape", NULL, NULL);
	if (window == NULL)
	{
		fprintf(stderr, "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials.\n");
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);

	// Initialize GLEW
	glewExperimental = true; // Needed in core profile
	if (glewInit() != GLEW_OK)
	{
		fprintf(stderr, "Failed to initialize GLEW\n");
		return -1;
	}

	// Ensure we can capture the escape key being pressed below
	glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);

	// Dark blue background
	glClearColor(0.0f, 0.0f, 0.4f, 0.0f);

	GLuint VertexArrayID;
	glGenVertexArrays(1, &VertexArrayID);
	glBindVertexArray(VertexArrayID);

	// Create and compile our GLSL program from the shaders
	GLuint programID = LoadShaders("vertexShader.vertexshader", "fragmentShader.fragmentshader");

	vector<GLfloat> vertexBufferData = vector<GLfloat>(dimX * dimX * 2, 0);
	unsigned int k = 0;
	for (unsigned int i = 0; i < dimX; ++i)
	{
		for (unsigned int j = 0; j < dimY; ++j)
		{
			vertexBufferData[k++] = (2 * (float)i / dimX); // X
			vertexBufferData[k++] = (2 * (float)j / dimY); // Y
		}
	}

	// This will identify our vertex buffer
	GLuint vertexbuffer;
	// Generate 1 buffer, put the resulting identifier in vertexbuffer
	glGenBuffers(1, &vertexbuffer);
	// The following commands will talk about our 'vertexbuffer' buffer
	glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
	// Give our vertices to OpenGL.
	glBufferData(GL_ARRAY_BUFFER, vertexBufferData.size(), vertexBufferData.begin()._Ptr, GL_STATIC_DRAW);

	unsigned int iteration = 0;

	// Generate weights
	vector<GLfloat> weightVector = vector<GLfloat>(dimX * dimY);

	do
	{
		unsigned int i = 0;
		for (const vector<TerrainTile>& v : field.getTerrain())
		{
			for (const TerrainTile& t : v)
			{
				weightVector[i++] = t.height;
			}
		}

		field.generateTerrain(1.5, .7, 2, iteration * .025, iteration * .2, 1, .2, .14);

		// Clear the screen
		glClear(GL_COLOR_BUFFER_BIT);

		// Use our shader
		glUseProgram(programID);

		// 1st attribute buffer : vertices
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
		glVertexAttribPointer(
			0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
			2,                  // size
			GL_FLOAT,           // type
			GL_FALSE,           // normalized?
			0,                  // stride
			(void*)0            // array buffer offset
		);

		// Draw
		for (unsigned int counter = 0; counter < dimX * dimY; ++counter)
		{
			glPointSize((GLfloat)10 * weightVector[counter]);
			glDrawArrays(GL_POINTS, counter, 1);
		}

		glDisableVertexAttribArray(0);

		// Swap buffers
		glfwSwapBuffers(window);
		glfwPollEvents();

		iteration += 1;
		Sleep(10);

	} // Check if the ESC key was pressed or the window was closed
	while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS &&
		glfwWindowShouldClose(window) == 0);

	// Cleanup VBO
	glDeleteBuffers(1, &vertexbuffer);
	glDeleteVertexArrays(1, &VertexArrayID);
	glDeleteProgram(programID);

	// Close OpenGL window and terminate GLFW
	glfwTerminate();

	return 0;
}