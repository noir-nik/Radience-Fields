#include <chrono>
#include <sstream>
#include "optimize.h"
#include "settings.h"
#include "save_bmp.h"

int enzyme_dup;
int enzyme_dupnoneed;
int enzyme_out;
int enzyme_const;

template <typename return_type, typename... T>
return_type __enzyme_autodiff(void *, T...);

template <typename return_type, typename... T>
return_type __enzyme_fwddiff(void *, T...);

void sh_eval_2(const float3 &rd, float *out)
{
	float x = rd.x, y = rd.y, z = rd.z, z2 = z * z;
	float c0, c1, s0, s1, tmp_a, tmp_b, tmp_c;

	out[0] = 0.28209479177387814;
	if (SH_WIDTH < 2) return;
	out[2] = z * 0.488602511902919923;
	out[6] = z2 * 0.94617469575756008 + -0.315391565252520045;
	c0 = x;
	s0 = y;

	tmp_a = -0.488602511902919978;
	out[3] = tmp_a * c0;
	out[1] = tmp_a * s0;
	tmp_b = z * -1.09254843059207896;
	out[7] = tmp_b * c0;
	out[5] = tmp_b * s0;
	c1 = x * c0 - y * s0;
	s1 = x * s0 + y * c0;

	tmp_c = 0.546274215296039478;
	out[8] = tmp_c * c1;
	out[4] = tmp_c * s1;
}



float eval_sh(float *sh, const float3& rayDir)
{
	float sh_coeffs[SH_WIDTH];
	sh_eval_2(rayDir, sh_coeffs);

	float sum = 0.0f;
	for (int i = 0; i < SH_WIDTH; i++)
		sum += sh[i] * sh_coeffs[i];

	return sum;
}

float2 RayBoxIntersection(float3 ray_pos, float3 ray_dir, float3 boxMin, float3 boxMax)
{
	ray_dir.x = 1.0f / ray_dir.x;
	ray_dir.y = 1.0f / ray_dir.y;
	ray_dir.z = 1.0f / ray_dir.z;

	float lo = ray_dir.x * (boxMin.x - ray_pos.x);
	float hi = ray_dir.x * (boxMax.x - ray_pos.x);

	float tmin = std::min(lo, hi);
	float tmax = std::max(lo, hi);

	float lo1 = ray_dir.y * (boxMin.y - ray_pos.y);
	float hi1 = ray_dir.y * (boxMax.y - ray_pos.y);

	tmin = std::max(tmin, std::min(lo1, hi1));
	tmax = std::min(tmax, std::max(lo1, hi1));

	float lo2 = ray_dir.z * (boxMin.z - ray_pos.z);
	float hi2 = ray_dir.z * (boxMax.z - ray_pos.z);

	tmin = std::max(tmin, std::min(lo2, hi2));
	tmax = std::min(tmax, std::max(lo2, hi2));

	return float2(tmin, tmax);
}

inline float RELU(float &val)
{
  return val > 0.0f ? val : 0.0f;
}

float3 renderVolumeTrilerp(Cell* grid, float tmin, float tmax, const float3 &rayPos, const float3 &rayDir, const RayMarcherExample* pImpl)
{

	float delta_t = 1.0f / pImpl->gridSize;
	float t = tmin;
	float log_light_intensity = 0.0f;

	float alpha = 1.0f;
	float3 color = float3(0.0f, 0.0f, 0.0f);

	while (t < tmax && alpha > 0.01f)
	{
		float3 position = rayPos + rayDir * t;

		float3 pos;
		pos.x = max(0.0f, min(position.x, static_cast<float>(pImpl->gridSize - 1)));
		pos.y = max(0.0f, min(position.y, static_cast<float>(pImpl->gridSize - 1)));
		pos.z = max(0.0f, min(position.z, static_cast<float>(pImpl->gridSize - 1)));

		int lx = max(0, min((int)(position.x * (float)pImpl->gridSize), pImpl->gridSize - 2));
		int ly = max(0, min((int)(position.y * (float)pImpl->gridSize), pImpl->gridSize - 2));
		int lz = max(0, min((int)(position.z * (float)pImpl->gridSize), pImpl->gridSize - 2));

		pos -= floor(pos);

		Cell &links000 = grid[lz * pImpl->gridSize * pImpl->gridSize + ly * pImpl->gridSize + lx];
		Cell &links001 = grid[lz * pImpl->gridSize * pImpl->gridSize + ly * pImpl->gridSize + (lx + 1)];
		Cell &links010 = grid[lz * pImpl->gridSize * pImpl->gridSize + (ly + 1) * pImpl->gridSize + lx];
		Cell &links011 = grid[lz * pImpl->gridSize * pImpl->gridSize + (ly + 1) * pImpl->gridSize + (lx + 1)];
		Cell &links100 = grid[(lz + 1) * pImpl->gridSize * pImpl->gridSize + ly * pImpl->gridSize + lx];
		Cell &links101 = grid[(lz + 1) * pImpl->gridSize * pImpl->gridSize + ly * pImpl->gridSize + (lx + 1)];
		Cell &links110 = grid[(lz + 1) * pImpl->gridSize * pImpl->gridSize + (ly + 1) * pImpl->gridSize + lx];
		Cell &links111 = grid[(lz + 1) * pImpl->gridSize * pImpl->gridSize + (ly + 1) * pImpl->gridSize + (lx + 1)];

		float3 rgb000, rgb001, rgb010, rgb011, rgb100, rgb101, rgb110, rgb111;
		rgb000.x = eval_sh(links000.sh_r, rayDir);
		rgb000.y = eval_sh(links000.sh_g, rayDir);
		rgb000.z = eval_sh(links000.sh_b, rayDir);
		rgb001.x = eval_sh(links001.sh_r, rayDir);
		rgb001.y = eval_sh(links001.sh_g, rayDir);
		rgb001.z = eval_sh(links001.sh_b, rayDir);
		rgb010.x = eval_sh(links010.sh_r, rayDir);
		rgb010.y = eval_sh(links010.sh_g, rayDir);
		rgb010.z = eval_sh(links010.sh_b, rayDir);
		rgb011.x = eval_sh(links011.sh_r, rayDir);
		rgb011.y = eval_sh(links011.sh_g, rayDir);
		rgb011.z = eval_sh(links011.sh_b, rayDir);
		rgb100.x = eval_sh(links100.sh_r, rayDir);
		rgb100.y = eval_sh(links100.sh_g, rayDir);
		rgb100.z = eval_sh(links100.sh_b, rayDir);
		rgb101.x = eval_sh(links101.sh_r, rayDir);
		rgb101.y = eval_sh(links101.sh_g, rayDir);
		rgb101.z = eval_sh(links101.sh_b, rayDir);
		rgb110.x = eval_sh(links110.sh_r, rayDir);
		rgb110.y = eval_sh(links110.sh_g, rayDir);
		rgb110.z = eval_sh(links110.sh_b, rayDir);
		rgb111.x = eval_sh(links111.sh_r, rayDir);
		rgb111.y = eval_sh(links111.sh_g, rayDir);
		rgb111.z = eval_sh(links111.sh_b, rayDir);


		float3& wb = pos;
		float3 wa = 1 - wb;

		float c00 = links000.density * wa.z + links001.density * wb.z;
		float c01 = links010.density * wa.z + links011.density * wb.z;
		float c10 = links100.density * wa.z + links101.density * wb.z;
		float c11 = links110.density * wa.z + links111.density * wb.z;
		float c0 = c00 * wa.y + c01 * wb.y;
		float c1 = c10 * wa.y + c11 * wb.y;
		float sigma = c0 * wa.x + c1 * wb.x;

		float3 col00, col01, col10, col11, col0, col1, rgb;
		col00 = rgb000 * wa.z + rgb001 * wb.z;
		col01 = rgb010 * wa.z + rgb011 * wb.z;
		col10 = rgb100 * wa.z + rgb101 * wb.z;
		col11 = rgb110 * wa.z + rgb111 * wb.z;
		col0 = col00 * wa.y + col01 * wb.y;
		col1 = col10 * wa.y + col11 * wb.y;
		rgb = col0 * wa.x + col1 * wb.x;

		float T = std::exp(log_light_intensity);

		float log_att = -RELU(sigma) * delta_t; 

		float weight = T * (1.0f - std::exp(log_att));

		color.x += weight * rgb.x;
		color.y += weight * rgb.y;
		color.z += weight * rgb.z;

		log_light_intensity += log_att;
		alpha = T;
		t += delta_t;
	}
	return color;
}

float3 renderVolume(Cell* grid, float tmin, float tmax, const float3 &rayPos, const float3 &rayDir, const RayMarcherExample* pImpl)
{

	float delta_t = 1.0f / pImpl->gridSize;
	float t = tmin;
	float log_light_intensity = 0.0f;

	float alpha = 1.0f;
	float3 color = float3(0.0f, 0.0f, 0.0f);

	while (t < tmax && alpha > 0.01f)
	{
		float3 position = rayPos + rayDir * t;

		int i = max(0, min(pImpl->gridSize - 1, (int)(position.x * (float)pImpl->gridSize))); 
		int j = max(0, min(pImpl->gridSize - 1, (int)(position.y * (float)pImpl->gridSize)));
		int k = max(0, min(pImpl->gridSize - 1, (int)(position.z * (float)pImpl->gridSize)));

		Cell &cell = grid[k * pImpl->gridSize * pImpl->gridSize + j * pImpl->gridSize + i];

		float sigma = cell.density;
		float T = std::exp(log_light_intensity);

		float log_att = -RELU(sigma) * delta_t;

		float weight = T * (1.0f - std::exp(log_att));

		float3 rgb(0.0f);
		rgb.x = eval_sh(cell.sh_r, rayDir);
		rgb.y = eval_sh(cell.sh_g, rayDir);
		rgb.z = eval_sh(cell.sh_b, rayDir);

		color += rgb * weight;

		log_light_intensity += log_att;
		alpha = T;
		t += delta_t;
	}
	return color;
}

void render(Cell* opt, float3* out, const RayMarcherExample* pImpl, const float4x4& c2w, const int W, const int H)
{
	float focal_inv = 1.0f / pImpl->focal;

# if USE_OMP
	#pragma omp parallel for
# endif
	for (int i = 0; i < H; i++){
		for (int j = 0; j < W; j++){
			float4 dir = {(j - W * 0.5f) * focal_inv,  (i - H * 0.5f) * focal_inv, -1.0f, 0.0f};
			float3 rayDir = to_float3(c2w * dir);
			float3 rayPos = to_float3(c2w.get_col(3));
			
			float2 tNearAndFar = RayBoxIntersection(rayPos, rayDir, pImpl->bb.min, pImpl->bb.max);

			out[i * W + j] = renderVolume(opt, tNearAndFar.x, tNearAndFar.y, rayPos, rayDir, pImpl);
		}
	}
}

float calculateMSE(Cell* opt, const float3* ref, const RayMarcherExample* pImpl, const float4x4& c2w)
{	
	float3* out = new float3[(pImpl->W) * (pImpl->H)];
	render(opt, out, pImpl, c2w, pImpl->W, pImpl->H);

	float resolution_inv = 1.0f / (pImpl->W * pImpl->H);
	float accum = 0.0;

	for (int i = 0; i < pImpl->W * pImpl->H; i++)
	{
		accum += ((out[i].x - ref[i].x) * (out[i].x - ref[i].x) + (out[i].y - ref[i].y) * (out[i].y - ref[i].y) + (out[i].z - ref[i].z) * (out[i].z - ref[i].z)) * resolution_inv;
	}
	delete[] out;
	return accum / NUM_REFFERENCES;
}

void RayMarcherExample::optimizeGrid(const float learning_rate, const int num_stages, const int* num_iterations, const float stop_threshold){
	std::chrono::high_resolution_clock::time_point start_stage, start_iter, end;
	std::chrono::duration<float> elapsed;
	Optimizer<float, OptimizerType::Adam> adam(learning_rate, momentum);
	bool converged = false;
	float mse;

	// Stages
	for (int i = 0; i < num_stages; i++){
		adam.setSize(gridSize * gridSize * gridSize * (sizeof(Cell) / sizeof(float)));
		std::cout << "Optimizing stage " << i + 1 << ", grid resolution: " << gridSize <<  std::endl;
		start_stage = std::chrono::high_resolution_clock::now();

		// Iterations
		for (int j = 0; j < num_iterations[i]; j++){
			printf("Iteration %d, ", j + 1);
			start_iter = std::chrono::high_resolution_clock::now();
			adam.zero_grad();
			for (auto& frame: frames){
				__enzyme_autodiff<void>((void*)calculateMSE,
					enzyme_dup, grid.data(), adam.ptr_grad(),
					enzyme_const, frame.image.data(),
					enzyme_const, this,
					enzyme_const, frame.transform_matrix);
				
			}

			// Step
			adam.update(grid.data());

			// Loss
			mse = calculateMSE(grid.data(), frames[0].image.data(), this, frames[0].transform_matrix);
			end = std::chrono::high_resolution_clock::now();
			elapsed = end - start_iter;
			std::cout << "MSE (image 0): " << mse << " iter_time: " << elapsed.count() << "s" << std::endl;
			if (mse < stop_threshold){
				std::cout << "Converged after " << j + 1 << " iterations, MSE: " << mse << std::endl;
				converged = true;
				break;
			}
			
			//Render
			int halfW = W / 2;
			int halfH = H / 2;
			std::vector<uint32_t> pixelData(halfW * halfH);
			std::vector<float3> realPixelData(halfW * halfH);
			focal *= 0.5f;
			render(grid.data(), realPixelData.data(), this, frames[1].transform_matrix, halfW, halfH);
			focal *= 2.0f;
			for (int ii = 0; ii < halfH; ii++)
				for (int jj = 0; jj < halfW; jj++)
					pixelData[ii * halfW + jj] = RealColorToUint32(to_float4(realPixelData[ii * halfW + jj], 1.0f));
			std::stringstream strOut;
			strOut << RENDER_OUTPUT_FOLDER << "stage_" << i + 1 << "_" << "iter_" << j + 1 << "_" << "size_" << gridSize << ".bmp";
			std::string fileName = strOut.str();
			generateBitmapImage(pixelData.data(), halfW, halfH, fileName.c_str(), true);
		}

		// End of stage
		end = std::chrono::high_resolution_clock::now();
		elapsed = end - start_stage;
		std::cout << "Stage " << i + 1 << " time: " << elapsed.count() << "s" << std::endl;
		
		if (converged) break;
		// Transition from smaller grid to larger grid
		if (i < num_stages - 1)
			subdivideGrid();
	}

}
