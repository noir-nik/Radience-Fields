#ifndef SETTINGS_H
#define SETTINGS_H

// ===== Folders =====
#define RENDER_OUTPUT_FOLDER "output/"
#define GRID_OUTPUT_FOLDER "output_grid_data/";

// ===== Use OpenMP =====
// If there is not enough RAM, disable it
#define USE_OMP 1

// ===== Number of references from dataset (max - 100) =====
#define NUM_REFFERENCES 20

// ===== Training settings =====
constexpr float learning_rate = 2.0f; // If Adam optimizer
//constexpr float learning_rate = 5.0e+5f; // If BasicGD optimizer

constexpr float momentum = 0.85f;
constexpr int num_stages = 2;
constexpr int num_iterations[num_stages] = {3, 1000};
constexpr int initialGridSize = 64;
constexpr int targetGridSize = 128;
constexpr float stop_threshold = 0.0005f;

#endif
