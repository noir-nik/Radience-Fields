#include "tracer.h"

#include <vector>
#include <chrono>
#include <string>

#include "tracer.h"

void sh_eval_2(const float3 &rd, float *out)
{
	float x = rd.x, y = rd.y, z = rd.z, z2 = z * z;
	float c0, c1, s0, s1, tmp_a, tmp_b, tmp_c;

	out[0] = 0.28209479177387814;
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

float eval_sh(float *sh, float3 rayDir)
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

static inline float3 EyeRayDir(float x, float y, float4x4 a_mViewProjInv)
{
	float4 pos = float4(2.0f * x - 1.0f, 2.0f * y - 1.0f, 0.0f, 1.0f);
	pos = a_mViewProjInv * pos;
	pos /= pos.w;
	return normalize(to_float3(pos));
}

static inline void transform_ray3f(float4x4 a_mWorldViewInv, float3 *ray_pos, float3 *ray_dir)
{
	float4 rayPosTransformed = a_mWorldViewInv * to_float4(*ray_pos, 1.0f);
	float4 rayDirTransformed = a_mWorldViewInv * to_float4(*ray_dir, 0.0f);

	(*ray_pos) = to_float3(rayPosTransformed);
	(*ray_dir) = to_float3(normalize(rayDirTransformed));
}

inline float RELU(float &val)
{
	return val > 0.0f ? val : 0.0f;
}



float4 RayMarcher::renderVolumeTrilerp(float tmin, float tmax, float3 &rayPos, float3 &rayDir)
{
	float delta_t = 1.0f / gridSize;
	float t = tmin;
	float log_light_intensity = 0.0f;

	float alpha = 1.0f;
	float4 color = float4(0.0f, 0.0f, 0.0f, 1.0f);

	while (t < tmax && alpha > 0.01f)
	{
		float3 position = rayPos + rayDir * t;

		float3 pos;
		pos.x = max(0.0f, min(position.x, static_cast<float>(gridSize - 1)));
		pos.y = max(0.0f, min(position.y, static_cast<float>(gridSize - 1)));
		pos.z = max(0.0f, min(position.z, static_cast<float>(gridSize - 1)));

		int lx = max(0, min((int)(position.x * (float)gridSize), gridSize - 2));
		int ly = max(0, min((int)(position.y * (float)gridSize), gridSize - 2));
		int lz = max(0, min((int)(position.z * (float)gridSize), gridSize - 2));
		pos -= floor(pos);

		Cell &links000 = grid[lz * gridSize * gridSize + ly * gridSize + lx];
		Cell &links001 = grid[lz * gridSize * gridSize + ly * gridSize + (lx + 1)];
		Cell &links010 = grid[lz * gridSize * gridSize + (ly + 1) * gridSize + lx];
		Cell &links011 = grid[lz * gridSize * gridSize + (ly + 1) * gridSize + (lx + 1)];
		Cell &links100 = grid[(lz + 1) * gridSize * gridSize + ly * gridSize + lx];
		Cell &links101 = grid[(lz + 1) * gridSize * gridSize + ly * gridSize + (lx + 1)];
		Cell &links110 = grid[(lz + 1) * gridSize * gridSize + (ly + 1) * gridSize + lx];
		Cell &links111 = grid[(lz + 1) * gridSize * gridSize + (ly + 1) * gridSize + (lx + 1)];

		float3 rgb000(eval_sh(links000.sh_r, rayDir), eval_sh(links000.sh_g, rayDir), eval_sh(links000.sh_b, rayDir));
		float3 rgb001(eval_sh(links001.sh_r, rayDir), eval_sh(links001.sh_g, rayDir), eval_sh(links001.sh_b, rayDir));
		float3 rgb010(eval_sh(links010.sh_r, rayDir), eval_sh(links010.sh_g, rayDir), eval_sh(links010.sh_b, rayDir));
		float3 rgb011(eval_sh(links011.sh_r, rayDir), eval_sh(links011.sh_g, rayDir), eval_sh(links011.sh_b, rayDir));
		float3 rgb100(eval_sh(links100.sh_r, rayDir), eval_sh(links100.sh_g, rayDir), eval_sh(links100.sh_b, rayDir));
		float3 rgb101(eval_sh(links101.sh_r, rayDir), eval_sh(links101.sh_g, rayDir), eval_sh(links101.sh_b, rayDir));
		float3 rgb110(eval_sh(links110.sh_r, rayDir), eval_sh(links110.sh_g, rayDir), eval_sh(links110.sh_b, rayDir));
		float3 rgb111(eval_sh(links111.sh_r, rayDir), eval_sh(links111.sh_g, rayDir), eval_sh(links111.sh_b, rayDir));

		float3 wb = pos;
		float3 wa = 1 - wb;

		float c00 = links000.density * wa.z + links001.density * wb.z;
		float c01 = links010.density * wa.z + links011.density * wb.z;
		float c10 = links100.density * wa.z + links101.density * wb.z;
		float c11 = links110.density * wa.z + links111.density * wb.z;
		float c0 = c00 * wa.y + c01 * wb.y;
		float c1 = c10 * wa.y + c11 * wb.y;
		float sigma = c0 * wa.x + c1 * wb.x;

		float3 col00 = rgb000 * wa.z + rgb001 * wb.z;
		float3 col01 = rgb010 * wa.z + rgb011 * wb.z;
		float3 col10 = rgb100 * wa.z + rgb101 * wb.z;
		float3 col11 = rgb110 * wa.z + rgb111 * wb.z;
		float3 col0 = col00 * wa.y + col01 * wb.y;
		float3 col1 = col10 * wa.y + col11 * wb.y;
		float3 rgb = col0 * wa.x + col1 * wb.x;


		float T = std::exp(log_light_intensity);

		float log_att = -RELU(sigma) * delta_t;

		float weight = T * (1.0f - std::exp(log_att));

		rgb = clamp(rgb, 0.0f, 1.0f);
		color += float4(rgb.x, rgb.y, rgb.z, 1.0f) * weight;


		log_light_intensity += log_att;
		alpha = T;
		t += delta_t;
	}
	return color;
}

static inline uint32_t RealColorToUint32(float4 real_color)
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

float4 RayMarcher::fragment(float2 fragCoord, uint2 iResolution)
{	
	
	float3 rayDir = EyeRayDir((float(fragCoord.x) + 0.5f) / float(iResolution.x), (float(fragCoord.y) + 0.5f) / float(iResolution.y), m_worldViewProjInv);
	float3 rayPos = float3(0.0f, 0.0f, 0.0f);
	transform_ray3f(m_worldViewInv, &rayPos, &rayDir);
	float2 tNearAndFar = RayBoxIntersection(rayPos, rayDir, bb.min, bb.max);

	float4 color = renderVolumeTrilerp(tNearAndFar.x, tNearAndFar.y, rayPos, rayDir);
	float4 fragColor = float4(pow(color.x, 1.0f / 2.2f), pow(color.y, 1.0f / 2.2f), pow(color.z, 1.0f / 2.2f), 1.0f);
	return fragColor;
}