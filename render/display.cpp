#include "settings.h"
#include "display.h"

#include <cstdio>

#define GLAD_GL_IMPLEMENTATION
#include <glad/glad.h> 
#include <GLFW/glfw3.h>

void error_callback(int error, const char* description)
{
	fprintf(stderr, "Error to: %s\n", description);
}

static GLFWwindow* window = NULL;

bool create_window(void){
	glfwSetErrorCallback(error_callback);
	if (!glfwInit()){
		fprintf(stderr, "Failed to initialize GLFW\n");
		return false;
	}
	
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_VISIBLE, GL_FALSE);

	window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "Ray-Marching", NULL, NULL);
	if (!window){
		fprintf(stderr, "Failed to create GLFW window\n");
		glfwTerminate();
		return false;
	}
	glfwMakeContextCurrent(window);
	
	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
		fprintf(stderr, "Failed to initialize GLAD\n");
		glfwTerminate();
		return false;
	}

	return true;
}

void destroy_window(void) {
	glfwDestroyWindow(window);
	glfwTerminate();
}
