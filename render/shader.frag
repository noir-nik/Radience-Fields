#version 450 core

#extension GL_GOOGLE_include_directive : enable
#include "settings.h"

uniform vec2 iResolution;
uniform mat4 uWorldViewProjInv;
uniform mat4 uWorldViewInv;
uniform int uGridSize;
uniform vec3 uBBoxMin;
uniform vec3 uBBoxMax;

out vec4 fragColor;

struct Cell {
	float density;
	float sh_r[SH_WIDTH];
	float sh_g[SH_WIDTH];
	float sh_b[SH_WIDTH];
};

layout(std430, binding = 3) buffer SSBO {
	Cell cells[];
};

void sh_eval_2(vec3 rd, out float sh_coeffs[9]) {
	float x = rd.x, y = rd.y, z = rd.z, z2 = z * z;
	float c0, c1, s0, s1, tmp_a, tmp_b, tmp_c;

	sh_coeffs[0] = 0.28209479177387814;
	sh_coeffs[2] = z * 0.488602511902919923;
	sh_coeffs[6] = z2 * 0.94617469575756008 + -0.315391565252520045;
	c0 = x;

	s0 = y;

	tmp_a = -0.488602511902919978;
	sh_coeffs[3] = tmp_a * c0;
	sh_coeffs[1] = tmp_a * s0;
	tmp_b = z * -1.09254843059207896;
	sh_coeffs[7] = tmp_b * c0;
	sh_coeffs[5] = tmp_b * s0;
	c1 = x * c0 - y * s0;
	s1 = x * s0 + y * c0;

	tmp_c = 0.546274215296039478;
	sh_coeffs[8] = tmp_c * c1;
	sh_coeffs[4] = tmp_c * s1;
}

float eval_sh(float sh[9], vec3 rayDir) {
	float sh_coeffs[9];
	sh_eval_2(rayDir, sh_coeffs);

	float sum = 0.0;
	for(int i = 0; i < 9; i++) sum += sh[i] * sh_coeffs[i];

	//return 1.0f * sh[0];
	return sum;
	//return clamp(sum, 0.0f, 1.0f);
}

vec2 RayBoxIntersection(vec3 ray_pos, vec3 ray_dir, vec3 boxMin, vec3 boxMax) {
	ray_dir.x = 1.0 / ray_dir.x;
	ray_dir.y = 1.0 / ray_dir.y;
	ray_dir.z = 1.0 / ray_dir.z;

	float lo = ray_dir.x * (boxMin.x - ray_pos.x);
	float hi = ray_dir.x * (boxMax.x - ray_pos.x);

	float tmin = min(lo, hi);

	float tmax = max(lo, hi);

	float lo1 = ray_dir.y * (boxMin.y - ray_pos.y);
	float hi1 = ray_dir.y * (boxMax.y - ray_pos.y);

	tmin = max(tmin, min(lo1, hi1));
	tmax = min(tmax, max(lo1, hi1));

	float lo2 = ray_dir.z * (boxMin.z - ray_pos.z);
	float hi2 = ray_dir.z * (boxMax.z - ray_pos.z);

	tmin = max(tmin, min(lo2, hi2));
	tmax = min(tmax, max(lo2, hi2));

	return vec2(tmin, tmax);
}

vec3 EyeRayDir(float x, float y) {
	vec4 pos = vec4(2.0 * x - 1.0, 2.0 * y - 1.0, 0.0, 1.0);
	pos = uWorldViewProjInv * pos;
	pos /= pos.w;
	return normalize(pos.xyz);
}

vec3 transform_ray3f(inout vec3 ray_pos, inout vec3 ray_dir) {
	vec4 rayPosTransformed = uWorldViewInv * vec4(ray_pos, 1.0);
	vec4 rayDirTransformed = uWorldViewInv * vec4(ray_dir, 0.0);

	ray_pos = vec3(rayPosTransformed);
	ray_dir = vec3(normalize(rayDirTransformed));

	return normalize(rayDirTransformed.xyz / rayDirTransformed.w);
}

float RELU(float val) {
	return max(val, 0.0);
}

vec4 renderVolumeTrilerp(float tmin, float tmax, vec3 rayPos, vec3 rayDir) {
    // Define delta_t
	float delta_t = 1.0 / uGridSize;
	float t = tmin;
	float log_light_intensity = 0.0;
	float alpha = 1.0;
	vec4 color = vec4(0.0, 0.0, 0.0, 1.0);

	while(t < tmax && alpha > 0.01) {
		vec3 position = rayPos + rayDir * t;
		vec3 pos = clamp(position, vec3(0.0), vec3(float(uGridSize - 1)));
		pos -= floor(pos);

		int lx = int(clamp(position.x * float(uGridSize), 0.0, float(uGridSize - 2)));
		int ly = int(clamp(position.y * float(uGridSize), 0.0, float(uGridSize - 2)));
		int lz = int(clamp(position.z * float(uGridSize), 0.0, float(uGridSize - 2)));

		Cell links000 = cells[lz * uGridSize * uGridSize + ly * uGridSize + lx];
		Cell links001 = cells[lz * uGridSize * uGridSize + ly * uGridSize + (lx + 1)];
		Cell links010 = cells[lz * uGridSize * uGridSize + (ly + 1) * uGridSize + lx];
		Cell links011 = cells[lz * uGridSize * uGridSize + (ly + 1) * uGridSize + (lx + 1)];
		Cell links100 = cells[(lz + 1) * uGridSize * uGridSize + ly * uGridSize + lx];
		Cell links101 = cells[(lz + 1) * uGridSize * uGridSize + ly * uGridSize + (lx + 1)];
		Cell links110 = cells[(lz + 1) * uGridSize * uGridSize + (ly + 1) * uGridSize + lx];
		Cell links111 = cells[(lz + 1) * uGridSize * uGridSize + (ly + 1) * uGridSize + (lx + 1)];

        // Calculate rgb values using spherical harmonics evaluation
		vec3 rgb000 = vec3(eval_sh(links000.sh_r, rayDir), eval_sh(links000.sh_g, rayDir), eval_sh(links000.sh_b, rayDir));
		vec3 rgb001 = vec3(eval_sh(links001.sh_r, rayDir), eval_sh(links001.sh_g, rayDir), eval_sh(links001.sh_b, rayDir));
		vec3 rgb010 = vec3(eval_sh(links010.sh_r, rayDir), eval_sh(links010.sh_g, rayDir), eval_sh(links010.sh_b, rayDir));
		vec3 rgb011 = vec3(eval_sh(links011.sh_r, rayDir), eval_sh(links011.sh_g, rayDir), eval_sh(links011.sh_b, rayDir));
		vec3 rgb100 = vec3(eval_sh(links100.sh_r, rayDir), eval_sh(links100.sh_g, rayDir), eval_sh(links100.sh_b, rayDir));
		vec3 rgb101 = vec3(eval_sh(links101.sh_r, rayDir), eval_sh(links101.sh_g, rayDir), eval_sh(links101.sh_b, rayDir));
		vec3 rgb110 = vec3(eval_sh(links110.sh_r, rayDir), eval_sh(links110.sh_g, rayDir), eval_sh(links110.sh_b, rayDir));
		vec3 rgb111 = vec3(eval_sh(links111.sh_r, rayDir), eval_sh(links111.sh_g, rayDir), eval_sh(links111.sh_b, rayDir));

		vec3 wa = 1 - pos;

        // Trilinear interpolation for density
		float c00 = links000.density * wa.z + links001.density * pos.z;
		float c01 = links010.density * wa.z + links011.density * pos.z;
		float c10 = links100.density * wa.z + links101.density * pos.z;
		float c11 = links110.density * wa.z + links111.density * pos.z;
		float c0 = c00 * wa.y + c01 * pos.y;
		float c1 = c10 * wa.y + c11 * pos.y;
		float sigma = c0 * wa.x + c1 * pos.x;

        // Trilinear interpolation for RGB
		vec3 col00 = rgb000 * wa.z + rgb001 * pos.z;
		vec3 col01 = rgb010 * wa.z + rgb011 * pos.z;
		vec3 col10 = rgb100 * wa.z + rgb101 * pos.z;
		vec3 col11 = rgb110 * wa.z + rgb111 * pos.z;
		vec3 col0 = col00 * wa.y + col01 * pos.y;
		vec3 col1 = col10 * wa.y + col11 * pos.y;
		vec3 rgb = col0 * wa.x + col1 * pos.x;

        // Compute attenuation
		float T = exp(log_light_intensity);
		float log_att = -RELU(sigma) * delta_t;
		float weight = T * (1.0 - exp(log_att));

        // Clamp
		rgb = clamp(rgb, 0.0, 1.0);
		color += vec4(rgb, 1.0) * weight;

		log_light_intensity += log_att;
		alpha = T;
		t += delta_t;
	}

	return color;
}

vec4 renderVolume(float tmin, float tmax, vec3 rayPos, vec3 rayDir) {

	float delta_t = 1.0 / uGridSize;
	float t = tmin;
	float log_light_intensity = 0.0;

	float alpha = 1.0;
	vec4 color = vec4(0.0, 0.0, 0.0, 1.0);

	while(t < tmax && alpha > 0.01) {
		vec3 position = rayPos + rayDir * t;

		int i = max(0, min(uGridSize - 1, int(position.x * float(uGridSize))));
		int j = max(0, min(uGridSize - 1, int(position.y * float(uGridSize))));
		int k = max(0, min(uGridSize - 1, int(position.z * float(uGridSize))));

		Cell cell = cells[k * uGridSize * uGridSize + j * uGridSize + i];

		float sigma = cell.density;

		float T = exp(log_light_intensity);

		float log_att = -RELU(sigma) * delta_t;

		float weight = T * (1.0 - exp(log_att));
		vec3 rgb = vec3(eval_sh(cell.sh_r, rayDir), eval_sh(cell.sh_g, rayDir), eval_sh(cell.sh_b, rayDir));
		
		// Clamp
		rgb = clamp(rgb, 0.0, 1.0);
		color += vec4(rgb, 0.0) * weight;

		log_light_intensity += log_att;
		alpha = T;
		t += delta_t;
	}
	return color;
}

void main() {

	vec3 rayDir = EyeRayDir((float(gl_FragCoord.x) + 0.5) / float(iResolution.x), (float(gl_FragCoord.y) + 0.5) / float(iResolution.y));
	vec3 rayPos = vec3(0.0, 0.0, 0.0);
	transform_ray3f(rayPos, rayDir);
	vec2 tNearAndFar = RayBoxIntersection(rayPos, rayDir, uBBoxMin, uBBoxMax);

	vec4 color = renderVolumeTrilerp(tNearAndFar.x, tNearAndFar.y, rayPos, rayDir);
	
	color = pow(color, vec4(1.0 / 2.2));
	fragColor = color;
}