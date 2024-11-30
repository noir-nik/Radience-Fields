#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <fstream>
#include <vector>
#include "LiteMath/LiteMath.h"
using namespace LiteMath;

struct shaderReader
{
	shaderReader(const std::string &filePath)
	{
		std::ifstream fileStream(filePath, std::ios::in | std::ios::binary);
		if (fileStream.is_open())
		{
			fileStream.seekg(0, std::ios::end);
			size_t size = fileStream.tellg();
			fileStream.seekg(0, std::ios::beg);

			source.resize(size + 1);
			fileStream.read(source.data(), size);
			source[size] = '\0';

			data = source.data();
			this->size = size;

			fileStream.close();
		}
		else
		{
			std::cerr << "Unable to open file: " << filePath << std::endl;
			exit(1);
		}
	}

	std::vector<char> source;
	const char *data;
	int size;
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

#endif
