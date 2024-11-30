#include "tracer.h"
#define STB_IMAGE_IMPLEMENTATION
#include <LiteMath/external/stb_image.h>

#include <json.hpp>
using json = nlohmann::json;

void RayMarcherExample::readJSONFromFile(const std::string& filename, size_t num_frames){
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Failed to open JSON file: " << filename << std::endl;
        exit(-1);
        return;
    }
    json j;
    file >> j;

    // Parse the additional integer and float values
    W = j["w"];
    H = j["h"];
    camera_angle_x = j["camera_angle_x"];

    focal = 0.5f * W / tan(0.5f * camera_angle_x);

    // Parse the JSON data into frames
    for (const auto& frame_data : j["frames"]) {
        if (num_frames > 0)
            if(frames.size() >= num_frames)
                break;
        
        Frame frame;
        frame.file_path = frame_data["file_path"];
        frame.rotation = frame_data["rotation"];

        // Parse the transform_matrix array into a float4x4 object
        auto transform_matrix_array = frame_data["transform_matrix"];
        float4x4 transform_matrix{
            transform_matrix_array[0][0], transform_matrix_array[0][1], transform_matrix_array[0][2], transform_matrix_array[0][3],
            transform_matrix_array[1][0], transform_matrix_array[1][1], transform_matrix_array[1][2], transform_matrix_array[1][3],
            transform_matrix_array[2][0], transform_matrix_array[2][1], transform_matrix_array[2][2], transform_matrix_array[2][3],
            transform_matrix_array[3][0], transform_matrix_array[3][1], transform_matrix_array[3][2], transform_matrix_array[3][3]
        };

        // XZY => XYZ // rows 1 <-> 2
        std::swap(transform_matrix.m_col[0][1], transform_matrix.m_col[0][2]);
        std::swap(transform_matrix.m_col[1][1], transform_matrix.m_col[1][2]);
        std::swap(transform_matrix.m_col[2][1], transform_matrix.m_col[2][2]);
        std::swap(transform_matrix.m_col[3][1], transform_matrix.m_col[3][2]);

        // Correct position (from Blender)
        transform_matrix = translate4x4(float3(0.5f, 0.33f, 0.5f)) * rotate4x4Y(180.0f * DEG_TO_RAD) * scale4x4(float3(0.4f, 0.4f, 0.4f)) * transform_matrix;

        frame.transform_matrix = transform_matrix;
        frames.push_back(frame);
    }
}

static void uint32_to_real_img(uint32_t* pixels, size_t w, size_t h, std::vector<float3>& img, bool flip_y = false){
	for (size_t y = 0; y < h; y++){
		for (size_t x = 0; x < w; x++){	
			img[y * h + x] = float3(pixels[y * h + x] & 0xFF, (pixels[y * h + x] >> 8) & 0xFF, (pixels[y * h + x] >> 16) & 0xFF) / 256.0f;
		}
	}
}

void RayMarcherExample::readImages(){
    int width, height, numChannels;
    stbi_set_flip_vertically_on_load(true);  
    for (auto& frame : frames) {
        unsigned char* imageData = stbi_load(frame.file_path.c_str(), &width, &height, &numChannels, 0);
        if (!imageData){
            std::cerr << "Error loading image file: " << frame.file_path.c_str() << std::endl;
            exit(-1);
            return;
        }
        if (width!= W || height != H) {
            std::cerr << "Error: image size mismatch" << std::endl;
            exit(-1);
            return;
        }
        
        // Uint32 to float3
        frame.image.resize(W * H);
        uint32_to_real_img((uint32_t*)imageData, W, H, frame.image);

        // free imageData
        stbi_image_free(imageData);
    }
}

void RayMarcherExample::subdivideGrid(size_t num_subdivisions){
    if (num_subdivisions == 0) return;
    std::vector<Cell> new_grid;
    std::cout << "Subdividing grid with " << num_subdivisions << " subdivisions from " << gridSize << " to " << gridSize * powl(2, num_subdivisions) << std::endl;
    
    for (size_t subdiv = 0; subdiv < num_subdivisions; subdiv++){

        new_grid.resize(gridSize * gridSize * gridSize * 8);

        for (size_t i = 0; i < gridSize; i++){
            for (size_t j = 0; j < gridSize; j++){
                for (size_t k = 0; k < gridSize; k++){
                    // 8 new cells
                    Cell& cell = grid[k * gridSize * gridSize + j * gridSize + i];
                    new_grid[2*k * gridSize * gridSize * 4 + 2*j* gridSize * 2 + 2*i] = cell;
                    new_grid[2*k * gridSize * gridSize * 4 + 2*j* gridSize * 2 + (2*i + 1)] = cell;
                    new_grid[2*k * gridSize * gridSize * 4 + (2*j+ 1) * gridSize * 2 + 2*i] = cell;
                    new_grid[2*k * gridSize * gridSize * 4 + (2*j+ 1) * gridSize * 2 + (2*i + 1)] = cell;
                    new_grid[(2*k + 1) * gridSize * gridSize * 4 + 2*j* gridSize * 2 + 2*i] = cell;
                    new_grid[(2*k + 1) * gridSize * gridSize * 4 + 2*j* gridSize * 2 + (2*i + 1)] = cell;
                    new_grid[(2*k + 1) * gridSize * gridSize * 4 + (2*j+ 1) * gridSize * 2 + 2*i] = cell;
                    new_grid[(2*k + 1) * gridSize * gridSize * 4 + (2*j+ 1) * gridSize * 2 + (2*i + 1)] = cell;
                }
            }
        }

        grid.resize(gridSize * gridSize * gridSize * 8);
        for (size_t i = 0; i < gridSize * gridSize * gridSize * 8; i++){
            grid[i] = new_grid[i];
        }
        gridSize *= 2;
    }    
}

void RayMarcherExample::loadGrid(const char* filePath, size_t _gridSize){
    
    gridSize = _gridSize;
    targetGridSize = _gridSize;
    grid.resize(gridSize * gridSize * gridSize);
    
    // Open the binary file for reading
    std::ifstream fin(filePath, std::ios::in | std::ios::binary);
    if (!fin.is_open()){
        std::cerr << "Cannot open file " << filePath << std::endl;
        exit(-1);
    }

    fin.read((char*)&grid[0], grid.size() * sizeof(Cell));
    if (!fin){
        std::cerr << "Error reading file " << filePath << std::endl;
        fin.close();
        exit(-1);
    }
    fin.close();
    std::cout << "Loaded grid from: " << filePath << std::endl;
}
void RayMarcherExample::saveGrid(const char *fileName, const size_t sh_width_target, size_t resolution){


    while (gridSize < resolution){
        std::cout << "Current grid size(" << gridSize << ") < resolution(" << resolution << ")" << std::endl;
        subdivideGrid();
    }
    if (gridSize > resolution){
        std::cout << "Current grid size(" << gridSize << ") > resolution(" << resolution << ")" << std::endl;
        std::cout << "Saving with resolution: " << gridSize << std::endl;
        resolution = gridSize;
    }

    std::vector<float> out_grid_data(resolution * resolution * resolution * (1 + 3 * sh_width_target), 0.0f);

    size_t min_width = SH_WIDTH < sh_width_target ? SH_WIDTH : sh_width_target;
    if(min_width != sh_width_target){
        std::cout << "Current SH_WIDTH(" << SH_WIDTH << ") != target SH_WIDTH(" << sh_width_target << ")" << std::endl;
        if (min_width < sh_width_target)
            std::cout << "SH " << min_width + 1 << " - " << sh_width_target << " are filled with 0" << std::endl;
        std::cout << "Saving with SH_WIDTH: " << sh_width_target << std::endl;
    }
    size_t out_cell_size = (1 + 3 * sh_width_target);
    for (size_t i = 0; i < resolution * resolution * resolution ; i++){
        out_grid_data[i * out_cell_size] = grid[i].density;
        for (size_t j = 0; j < min_width; j++){
            out_grid_data[i * out_cell_size + 1 + j] = grid[i].sh_r[j];
            out_grid_data[i * out_cell_size + 1 + sh_width_target + j] = grid[i].sh_g[j];
            out_grid_data[i * out_cell_size + 1 + 2 * sh_width_target + j] = grid[i].sh_b[j];
        }
    }
	
	
	std::ofstream outfile(fileName, std::ios::binary);
	outfile.write(reinterpret_cast<char*>(out_grid_data.data()), out_grid_data.size() * sizeof(float));
	outfile.close();
    
}
