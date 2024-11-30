#include "save_bmp.h"

#define USE_ALPHA 0

#if USE_ALPHA
const int BYTES_PER_PIXEL = 4; // Red, green, blue, & alpha
#else
const int BYTES_PER_PIXEL = 3; // Red, green, blue
#endif
const int FILE_HEADER_SIZE = 14;
const int INFO_HEADER_SIZE = 40;

void generateBitmapImage(const uint32_t* image, int height, int width, const char* imageFileName, bool flip) {
    int widthInBytes = width * BYTES_PER_PIXEL;

    unsigned char padding[3] = {0, 0, 0};
    int paddingSize = (4 - (widthInBytes) % 4) % 4;

    int stride = (widthInBytes) + paddingSize;

    FILE* imageFile = fopen(imageFileName, "wb");
    if (!imageFile) {
        fprintf(stderr, "Error: Unable to open file %s for writing\n", imageFileName);
        return;
    }

    unsigned char* fileHeader = createBitmapFileHeader(height, stride);
    fwrite(fileHeader, 1, FILE_HEADER_SIZE, imageFile);

    unsigned char* infoHeader = createBitmapInfoHeader(height, width);
    fwrite(infoHeader, 1, INFO_HEADER_SIZE, imageFile);

	struct Pixel { 
        unsigned char r, g, b;
        #if USE_ALPHA
            unsigned char a; 
        #endif
    };
	std::vector<Pixel> pixels2(width * height);
	for (size_t i = 0; i < pixels2.size(); i++)
	{
		Pixel px;
		px.r       = (image[i] & 0x00FF0000) >> 16;
		px.g       = (image[i] & 0x0000FF00) >> 8;
		px.b       = (image[i] & 0x000000FF);
        #if USE_ALPHA
		    px.a       = (image[i] & 0xFF000000) >> 24;
        #endif
		pixels2[i] = px;
	}

	if (flip) {
		for (int i = 0; i < height; ++i) {
			fwrite(&pixels2[0] + (i * width), BYTES_PER_PIXEL, width, imageFile);
			fwrite(padding, 1, paddingSize, imageFile);
    	}
	} else {
		for (int i = height - 1; i >= 0; --i) {
			fwrite(&pixels2[0] + (i * width), BYTES_PER_PIXEL, width, imageFile);
			fwrite(padding, 1, paddingSize, imageFile);
		}
	}

    fclose(imageFile);
}

unsigned char* createBitmapFileHeader(int height, int stride) {
    int fileSize = FILE_HEADER_SIZE + INFO_HEADER_SIZE + (stride * height);

    static unsigned char fileHeader[] = {
        0, 0,     // Signature 'BM'
        0, 0, 0, 0, // Image file size in bytes
        0, 0, 0, 0, // Reserved
        0, 0, 0, 0, // Start of pixel array
    };

    fileHeader[0] = 'B';
    fileHeader[1] = 'M';
    fileHeader[2] = (unsigned char)(fileSize);
    fileHeader[3] = (unsigned char)(fileSize >> 8);
    fileHeader[4] = (unsigned char)(fileSize >> 16);
    fileHeader[5] = (unsigned char)(fileSize >> 24);
    fileHeader[10] = FILE_HEADER_SIZE + INFO_HEADER_SIZE;

    return fileHeader;
}

unsigned char* createBitmapInfoHeader(int height, int width) {
    static unsigned char infoHeader[] = {
        0, 0, 0, 0, // Header size
        0, 0, 0, 0, // Image width
        0, 0, 0, 0, // Image height
        0, 0,       // Number of color planes
        0, 0,       // Bits per pixel
        0, 0, 0, 0, // Compression
        0, 0, 0, 0, // Image size
        0, 0, 0, 0, // Horizontal resolution
        0, 0, 0, 0, // Vertical resolution
        0, 0, 0, 0, // Colors in color table
        0, 0, 0, 0, // Important color count
    };

    infoHeader[0] = INFO_HEADER_SIZE;
    infoHeader[4] = width;
    infoHeader[5] = width >> 8;
    infoHeader[6] = width >> 16;
    infoHeader[7] = width >> 24;
    infoHeader[8] = height;
    infoHeader[9] = height >> 8;
    infoHeader[10] = height >> 16;
    infoHeader[11] = height >> 24;
    infoHeader[12] = 1;  // Number of color planes
    infoHeader[14] = BYTES_PER_PIXEL * 8; // Bits per pixel (32 for RGBA)

    return infoHeader;
}