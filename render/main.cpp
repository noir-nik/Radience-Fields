#include "settings.h"
#include "display.h"
#include "tracer.h"
#include "utils.hpp"
#include <csignal>
#include <chrono>
#include <sstream>
#include <iomanip>
#include <glad/glad.h>
#include <memory>

#include"save_bmp.h"

#define OUTPUT_FOLDER "output/"

bool is_running = true;

constexpr int SCREEN_WIDTH = WINDOW_WIDTH;
constexpr int SCREEN_HEIGHT = WINDOW_HEIGHT;

int main(int argc, char **argv)
{	
	signal(SIGINT, [](int){exit(0);});
	std::chrono::high_resolution_clock::time_point start, start_copy, end;
	std::chrono::duration<float> elapsed;
	std::shared_ptr<RayMarcher> pImpl = nullptr;
	std::vector<uint32_t> pixelData(SCREEN_WIDTH * SCREEN_HEIGHT);

	if (argc == 1 || argc > 5 || (argc > 1 && strcmp(argv[1], "--help") == 0))
	{
		printf("Usage: %s [--help] [--m <path to binary grid> ] [--s <grid size> ]\n", "rf-render");
		return 0;
	}

	int gridSize = 128;
	std::string dataFilePath = "data/model.dat";
	for (int i = 1; i < argc; ++i)
	{
		if (strcmp(argv[i], "--model") == 0 || strcmp(argv[i], "-m") == 0)
		{
			dataFilePath = argv[++i];
		}
		else if (strcmp(argv[i], "--size") == 0 || strcmp(argv[i], "-s") == 0)
		{
			gridSize = atoi(argv[++i]);
		}
	}

	pImpl = std::make_shared<RayMarcher>();
	pImpl->InitGrid(gridSize);
	pImpl->SetBoundingBox(float3(0, 0, 0), float3(1, 1, 1));
	pImpl->SetWorldViewMProjatrix(perspectiveMatrix(45, 1, 0.1, 100));

	// Read binary grid
	std::ifstream fin(dataFilePath, std::ios::in | std::ios::binary);
	if (!fin.is_open()){
		std::cerr << "Cannot open file " << dataFilePath << std::endl;
		return -1;
	}
	fin.read((char*)&pImpl->grid[0], pImpl->grid.size() * sizeof(Cell));
	if (!fin){
		std::cerr << "Error reading file " << dataFilePath << std::endl;
		fin.close();
		return -1;
	}
	fin.close();

	/*
	//////////////// CPU //////////////////
	*/
	printf("CPU render starts...\n");
	for(int k = 0; k < 7; k++)
	{	
		float3 pos = float3(0.0, 0.0, 1.3);
		float4x4 viewMat = lookAt(pos, float3(0.0, 0.0, 0.0), float3(0.0, 1.0, 0.0)) * rotate4x4Y(-float(360.0 / 7 * k) * DEG_TO_RAD) * translate4x4(float3(-0.5, -0.5, -0.5));
		pImpl->SetWorldViewMatrix(viewMat);
		start = std::chrono::high_resolution_clock::now();
#pragma omp parallel for
		for (int i = 0; i < SCREEN_HEIGHT; i++){
			for (int j = 0; j < SCREEN_WIDTH; j++){
				float4 color_vec = pImpl->fragment(float2 (j, i), uint2(SCREEN_WIDTH, SCREEN_HEIGHT));
				pixelData[i * SCREEN_WIDTH + j] = RealColorToUint32(color_vec);
			}
		}
		end = std::chrono::high_resolution_clock::now();
		elapsed = end-start;
		std::cout << "img no. = " << k + 1 << ", timeRender = " << elapsed.count() * 1000 << " ms" << std::endl;

		std::stringstream strOut;
		strOut << std::fixed << std::setprecision(2) << OUTPUT_FOLDER << "out_cpu_" << k + 1 << ".bmp";
		std::string fileName = strOut.str();

		generateBitmapImage(pixelData.data(), SCREEN_WIDTH, SCREEN_HEIGHT, fileName.c_str(), true);
	}
	/*
	//////////////// GPU //////////////////
	*/
	int success;
	char infoLog[512];
	if (!create_window())
	{
		std::cerr << "Failed to create window" << std::endl;
		return -1;
	}
	// Vertex shader
	const char *vertexShaderSource = R"(
		#version 450
		layout (location = 0) in vec2 aPos;
		void main() {
			gl_Position = vec4(aPos.x, aPos.y, 0.0, 1.0);
		}
	)";

	GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
	glCompileShader(vertexShader);

	glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
	if (!success) {
		glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
		std::cerr << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n" << infoLog << std::endl;
		return -1;
	}

	// Fragment shader
	GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
	shaderReader shaderFile(FRAGMENT_SOURCE);
	glShaderSource(fragmentShader, 1, &shaderFile.data, NULL);
	glCompileShader(fragmentShader);
	glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
	if (!success) {
		glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog);
		std::cerr << "ERROR::SHADER::FRAGMENT::COMPILATION_FAILED\n" << infoLog << std::endl;
		return -1;
	}

	// Shader program
	GLuint shaderProgram = glCreateProgram();
	glAttachShader(shaderProgram, vertexShader);
	glAttachShader(shaderProgram, fragmentShader);
	glLinkProgram(shaderProgram);
	glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
	if (!success) {
		glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
		std::cerr << "ERROR::SHADER::PROGRAM::LINKING_FAILED\n" << infoLog << std::endl;
		return -1;
	}
	glDetachShader(shaderProgram, vertexShader);
	glDetachShader(shaderProgram, fragmentShader);
	glDeleteShader(vertexShader);
	glDeleteShader(fragmentShader);

	glUseProgram(shaderProgram);
	printf("GPU render starts...\n");
	glFinish();
	start_copy = std::chrono::high_resolution_clock::now();

		
	// Uniforms
	GLint resolutionLoc = glGetUniformLocation(shaderProgram, "iResolution");
	GLint worldViewProjInvLoc = glGetUniformLocation(shaderProgram, "uWorldViewProjInv");
	GLint worldViewInvLoc = glGetUniformLocation(shaderProgram, "uWorldViewInv");
	GLint uGridSizeLoc = glGetUniformLocation(shaderProgram, "uGridSize");
	GLint uBBoxMinLoc = glGetUniformLocation(shaderProgram, "uBBoxMin");
	GLint uBBoxMaxLoc = glGetUniformLocation(shaderProgram, "uBBoxMax");

	// Vertex data
	GLfloat vertices[] = {
		-1.0f,  1.0f,
		-1.0f, -1.0f,
		1.0f, -1.0f,
		1.0f,  1.0f
	};
	
	unsigned int indices[] = {
		0, 1, 2,
		0, 2, 3
	};

	// VAO
	GLuint VAO, VBO, EBO;
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);
	glGenBuffers(1, &EBO);

	glBindVertexArray(VAO);

	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);

	glBindVertexArray(0);

	///// Copy data to GPU /////
	glUseProgram(shaderProgram);
	glBindVertexArray(VAO);

	// Grid data
	GLuint gridBuffer;
	glGenBuffers(1, &gridBuffer);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, gridBuffer);
	glBufferData(GL_SHADER_STORAGE_BUFFER, gridSize * gridSize * gridSize * sizeof(Cell), pImpl->grid.data(), GL_STATIC_DRAW);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, gridBuffer);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
	
	// Uniforms
	glUniform3f(uBBoxMinLoc, pImpl->bb.min.x, pImpl->bb.min.y, pImpl->bb.min.z);
	glUniform3f(uBBoxMaxLoc, pImpl->bb.max.x, pImpl->bb.max.y, pImpl->bb.max.z);
	glUniform1i(uGridSizeLoc, gridSize);
	glUniform2f(resolutionLoc, SCREEN_WIDTH, SCREEN_HEIGHT);
	glUniformMatrix4fv(worldViewProjInvLoc, 1, GL_FALSE, &(pImpl->m_worldViewProjInv[0][0]));

	glFinish();
	start = std::chrono::high_resolution_clock::now();
	end = std::chrono::high_resolution_clock::now();
	elapsed = end - start;
	std::cout << "GPU copy time: " << elapsed.count()  * 1000 << "ms" << std::endl;

	for (int k = 0; k < 7; k++) {
		float3 pos = float3(0.0, 0.0, 1.3);
		float4x4 viewMat = lookAt(pos, float3(0.0, 0.0, 0.0), float3(0.0, 1.0, 0.0)) * rotate4x4Y(-float(360.0 / 7 * k) * DEG_TO_RAD) * translate4x4(float3(-0.5, -0.5, -0.5));
		pImpl->SetWorldViewMatrix(viewMat);
		glUniformMatrix4fv(worldViewInvLoc, 1, GL_FALSE, &(pImpl->m_worldViewInv[0][0]));
		glFinish();
		start = std::chrono::high_resolution_clock::now();
		glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
		glFinish();
		end = std::chrono::high_resolution_clock::now();
		elapsed = end - start;
		std::cout << "img no. = " << k + 1 << ", GPU render time: "  << elapsed.count() * 1000 << "ms" << std::endl;
		glReadBuffer(GL_BACK);
		glPixelStorei(GL_PACK_ALIGNMENT, 1);
		glReadPixels(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT, GL_RGBA, GL_UNSIGNED_BYTE, pixelData.data());
		std::stringstream strOut;
		strOut << std::fixed << std::setprecision(2) << OUTPUT_FOLDER << "out_gpu_" << k + 1 << ".bmp";
		std::string fileName = strOut.str();
		generateBitmapImage(pixelData.data(), SCREEN_WIDTH, SCREEN_HEIGHT, fileName.c_str(), true);
	}	

	// Cleanup
	glDeleteVertexArrays(1, &VAO);
	glDeleteBuffers(1, &VBO);
	glDeleteBuffers(1, &EBO);
	glDeleteBuffers(1, &gridBuffer);
	glDeleteProgram(shaderProgram);

	destroy_window();
	return 0;
}