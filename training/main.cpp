#include <csignal>
#include <unordered_set>

#include "settings.h"
#include "includes.h"
#include "save_bmp.h"
#include "optimize.h"

const std::unordered_set<std::string> training_sets = {"lego", "chair"};
std::shared_ptr<RayMarcherExample> pImpl = nullptr;
bool initialized = false;
int main(int argc, char **argv)
{	
	signal(SIGINT, [](int){
		if (initialized){
			std::cout << "Save grid? (y)" << std::endl;
			std::string save;
			std::cin >> save;
			if (save == "y" || save == "Y"){
				std::cout << "Saving grid..." << std::endl;
				std::stringstream strOut;
				strOut << GRID_OUTPUT_FOLDER"tmp_grid_" << pImpl->gridSize << "_sh_" << SH_WIDTH << ".bin";
				std::string fileName = strOut.str();
				pImpl->saveGrid(fileName.c_str(), 9, pImpl->gridSize);
				//pImpl->saveGrid(fileName.c_str(), SH_WIDTH, pImpl->gridSize);
			}
		}
		exit(0);
	});
	std::chrono::high_resolution_clock::time_point start, start_copy, end;
	std::chrono::duration<float> elapsed;
	
	pImpl = std::make_shared<RayMarcherExample>();
	std::vector<uint32_t> pixelData;

    /*
	///////////////// Choose training set //////////////////
	*/
	std::cout << "Available training sets:" << std::endl;
	for (const auto& training_set : training_sets) {
		std::cout << training_set << std::endl;
	}
    std::string training_set_name;
	while (training_set_name.empty()){
		std::cout << "Enter the name of the training set: ";
		std::cin >> training_set_name;
		if (training_set_name.size() <= 2) exit(0);
		if (training_sets.find(training_set_name) == training_sets.end()) {
			std::cerr << "Error: Invalid training set name." << std::endl;
			training_set_name.clear();
		}
	}
	std::string transforms_file = "train/" + training_set_name + "/transforms.json";

	/*
	///////////////// Read JSON //////////////////
	*/
	pImpl->readJSONFromFile(transforms_file, NUM_REFFERENCES);
	std::cout << "Read JSON from file: " << transforms_file << std::endl;
	pImpl->readImages();
	std::cout << "Read images: " << pImpl->frames[0].file_path << " ~ " << pImpl->frames.back().file_path << std::endl;
	const int& SCREEN_WIDTH = pImpl->W;
	const int& SCREEN_HEIGHT = pImpl->H;
	pixelData.resize(SCREEN_WIDTH * SCREEN_HEIGHT);

	/*
	/////////////// Init grid //////////////////
	*/
	
	pImpl->InitGrid(initialGridSize, targetGridSize);
	pImpl->SetBoundingBox(float3(0, 0, 0), float3(1, 1, 1));
	pImpl->SetWorldViewMProjatrix(perspectiveMatrix(45, 1, 0.1, 100));
	initialized = true;

	/*
	/////////////// Optimize grid //////////////////
	*/
	start = std::chrono::high_resolution_clock::now();
	pImpl->optimizeGrid(learning_rate, num_stages, num_iterations, stop_threshold);
	end = std::chrono::high_resolution_clock::now();
	elapsed = end - start;
	std::cout << "Optimize time: " << elapsed.count() << std::endl;

	/*
	/////////////// Render //////////////////
	*/
	for(int k = 0; k < pImpl->frames.size(); k++){	
		start = std::chrono::high_resolution_clock::now();
		std::vector<float3> realPixelData(SCREEN_WIDTH * SCREEN_HEIGHT);
		render(pImpl->grid.data(), realPixelData.data(), pImpl.get(), pImpl->frames[k].transform_matrix, pImpl->W, pImpl->H);
		for (int i = 0; i < SCREEN_HEIGHT; i++)
			for (int j = 0; j < SCREEN_WIDTH; j++)
				pixelData[i * SCREEN_WIDTH + j] = RealColorToUint32(to_float4(realPixelData[i * SCREEN_WIDTH + j], 1.0f));
		end = std::chrono::high_resolution_clock::now();
		elapsed = end-start;
		std::cout << "img no. = " << k << ", timeRender = " << elapsed.count() * 1000 << " ms" << std::endl;
		std::stringstream strOut;
		strOut << RENDER_OUTPUT_FOLDER << "out_" << training_set_name << "_" << k << ".bmp";
		std::string fileName = strOut.str();
		generateBitmapImage(pixelData.data(), SCREEN_WIDTH, SCREEN_HEIGHT, fileName.c_str(), true);
	}
	
	/*
	/////////////// Save grid //////////////////
	*/
	// if output folder doesnt exist create it
	if (mkdir(GRID_OUTPUT_FOLDER, 0777) == -1) {
		if (errno != EEXIST) {
			std::cerr << "Error creating folder: " << GRID_OUTPUT_FOLDER << std::endl;
			return 1;
		}
	}
	std::stringstream strOut;
	strOut << GRID_OUTPUT_FOLDER<< training_set_name << "_grid_" << pImpl->gridSize << ".bin";
	std::string fileName = strOut.str();
	pImpl->saveGrid(fileName.c_str(), 9, 128);

	return 0;
}



