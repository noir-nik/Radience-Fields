#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <cstdint>

#include "LiteMath/LiteMath.h"
using namespace LiteMath;

const size_t SH_WIDTH = 1;

struct Frame {
    std::string file_path;
    float rotation;
    float4x4 transform_matrix;
	std::vector<float3> image;
};

struct Cell
{
	float density;
	float sh_r[SH_WIDTH];
	float sh_g[SH_WIDTH];
	float sh_b[SH_WIDTH];

	// operator =
	Cell& operator=(const Cell& other)
	{
		density = other.density;
		for (size_t i = 0; i < SH_WIDTH; i++)
		{
			sh_r[i] = other.sh_r[i];
			sh_g[i] = other.sh_g[i];
			sh_b[i] = other.sh_b[i];
		}
		return *this;
	}
};

struct BoundingBox
{
	float3 min;
	float3 max;
};

class RayMarcherExample // : public IRenderAPI
{
public:
	RayMarcherExample()
	{
		const float4x4 view = lookAt(float3(0, 1.5, -3), float3(0, 0, 0), float3(0, 1, 0)); // pos, look_at, up
		const float4x4 proj = perspectiveMatrix(90.0f, 1.0f, 0.1f, 100.0f);
		m_worldViewInv = inverse4x4(view);
		m_worldViewProjInv = inverse4x4(proj);
	}

	void InitGrid(const int _initialGridSize, const float _targetGridSize)
	{	
		gridSize = _initialGridSize;
		targetGridSize = _targetGridSize;
		grid.resize(gridSize * gridSize * gridSize);

		for (size_t i = 0; i < gridSize * gridSize * gridSize; i++)
		{
			grid[i].density = 0.01f;
			for (size_t j = 0; j < SH_WIDTH; j++)
			{
				grid[i].sh_r[j] = 0.1f;
				grid[i].sh_g[j] = 0.1f;
				grid[i].sh_b[j] = 0.1f;
			}
		}

	}
	void loadGrid(const char* filePath, size_t _gridSize);
	void saveGrid(const char *fileName, const size_t sh_width = 9, size_t resolution = 128);

	void subdivideGrid(size_t num_subdivisions = 1);

	void SetBoundingBox(const float3 boxMin, const float3 boxMax)
	{
		bb.min = boxMin;
		bb.max = boxMax;
	}

	void SetWorldViewMProjatrix(const float4x4 &a_mat) { m_worldViewProjInv = inverse4x4(a_mat); }
	void SetWorldViewMatrix(const float4x4 &a_mat) { m_worldViewInv = inverse4x4(a_mat); }

	float4 fragment(float2 fragCoord, uint2 iResolution);
	float4 renderVolume(float tmin, float tmax, float3& rayPos, float3& rayDir);
	float4 renderVolumeTrilerp(float tmin, float tmax, float3 &rayPos, float3 &rayDir);
	
	void readJSONFromFile(const std::string& filename, size_t num_frames = 0);
	void readImages();

	void optimizeGrid(const float learning_rate, const int num_stages, const int* num_iterations, const float stop_threshold);

	std::vector<Cell> grid;
	int gridSize;
	int targetGridSize;
	BoundingBox bb;
	
	float rayMarchTime;

	std::vector<Frame> frames;
	float4x4 m_worldViewProjInv;
	float4x4 m_worldViewInv;
	int W, H;
    float camera_angle_x;
	float focal;
};

inline uint32_t RealColorToUint32(float4 real_color)
{
	float r = real_color[0] * 256.0f;
	float g = real_color[1] * 256.0f;
	float b = real_color[2] * 256.0f;
	float a = real_color[3] * 256.0f;

	uint32_t red = max(0, min(255, (int)r));
	uint32_t green = max(0, min(255, (int)g));
	uint32_t blue = max(0, min(255, (int)b));
	uint32_t alpha = max(0, min(255, (int)a));

	return red | (green << 8) | (blue << 16) | (alpha << 24);
}
