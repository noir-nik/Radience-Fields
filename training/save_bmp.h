#ifndef SAVE_BMP
#define SAVE_BMP
#include <fstream>
#include <iostream>
#include <cstdint>
#include <cstring>
#include <vector>

unsigned char* createBitmapInfoHeader(int height, int width);
unsigned char* createBitmapFileHeader(int height, int stride);
void generateBitmapImage(const uint32_t* image, int height, int width, const char* imageFileName, bool flip = false);

#endif