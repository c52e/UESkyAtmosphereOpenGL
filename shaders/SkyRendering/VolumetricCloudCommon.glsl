#ifndef _VOLUMETRIC_CLOUD_COMMON_GLSL
#define _VOLUMETRIC_CLOUD_COMMON_GLSL

#include "AtmosphereInterface.glsl"

layout(std140, binding = 1) uniform VolumetricCloudCommonBufferData{
	mat4 uInvMVP;
	mat4 uReprojectMat;
	mat4 uReprojectMatNoProjection;

	mat4 uLightVP;
	mat4 uInvLightVP;
	mat4 uShadowMapReprojectMat;
	mat4 uScreenToShadowMap;

	vec3 uCameraPos;
	uint uBaseShadingIndex;

	float uZNear;
	float uLightFroxelMaxDist;
	float uBottomAltitude;
	float uTopAltitude;

	vec3 uSunDirection;
	float uFrameID;

	float uInvShadowFroxelMaxDistance;
	float uAerialPerspectiveLutMaxDistance;
	float uShadowFroxelMaxDistance;
	float uEarthRadius;

};

const float kMinTransmittance = 0.01;

float DepthToLinearDepth(float depth) {
	return min(uZNear / depth, 1e9);
}

float CalHeight01(vec3 pos) {
	float altitude = length(vec3(pos.xy, pos.z + uEarthRadius)) - uEarthRadius;
	return (altitude - uBottomAltitude) / (uTopAltitude - uBottomAltitude);
}

// Same order as textureGather
ivec2 IndexToOffset(uint index) {
#if 0
	const ivec2 offsets[4] = {{0, 1}, {1, 1}, {1, 0}, {0, 0}};
	return offsets[index];
#else
	return ivec2(
		((index + 1) >> 1) & 1,
		((index + 2) >> 1) & 1
	);
#endif
}

vec4 texelFetchClamp(sampler2D s, ivec2 p, int lod) {
	return texelFetch(s, clamp(p, ivec2(0), textureSize(s, lod) - 1), lod);
}

float HenyeyGreensteinInvertcdf(float xi, float g) {
	float one_plus_g2 = 1.0 + g * g;
	float one_minus_g2 = 1.0 - g * g;
	float one_over_2g = 0.5 / g;
	float t = (one_minus_g2) / (1.0 - g + 2.0 * g * xi);
	return one_over_2g * (one_plus_g2 - t * t);
}

vec3 GetSunVisibility(sampler2D transmittance_texture, vec3 pos) {
	vec3 up_dir = vec3(pos.xy, pos.z + uEarthRadius);
	float r = length(up_dir);
	up_dir /= r;
	float mu_s = dot(uSunDirection, up_dir);
	return GetSunVisibility(transmittance_texture, r, mu_s);
}

vec3 GetAerialPerspective(sampler3D aerial_perspective_luminance_texture,
		sampler3D aerial_perspective_transmittance_texture,
		vec2 uv, float t, float r, float mu, out vec3 transmittance) {
	if (r > top_radius) {
		float near_distance;
		if (FromSpaceIntersectTopAtmosphereBoundary(r, mu, near_distance)) {
			t -= near_distance;
		}
		else {
			t = 0;
		}
	}
	vec3 uvw = GetAerialPerspectiveTextureUvwFromTexCoordDistance(
		uv, t, uAerialPerspectiveLutMaxDistance, textureSize(aerial_perspective_luminance_texture, 0));
	transmittance = texture(aerial_perspective_transmittance_texture, uvw).rgb;
	vec3 luminance = texture(aerial_perspective_luminance_texture, uvw).rgb;
	return luminance;
}

#define CLOUD_SAMPLE_MAIN_RAY 0
#define CLOUD_SAMPLE_SHADOW_RAY 1
#define CLOUD_SAMPLE_SHADOW_MAP 2
#define CLOUD_SAMPLE_PATH_TRACING 3

struct CloudSampleParams {
	vec3 pos;
	float height01;
	vec3 camera_pos;
	int caller_type;
};

struct CloudSampleResult {
	float sigma_t;
	float phase_g;

	// Oz: The Great and Volumetric - Magnus Wrenninge
	float attenuation;
	float contribution;
};

CloudSampleResult DefaultCloudSampleResult() {
	CloudSampleResult res;
	res.sigma_t = 0;
	res.phase_g = 0.85;
	res.attenuation = 0.5;
	res.contribution = 0.5;
	return res;
}

// Should be defined in material shader (if used)
void SampleCloud(CloudSampleParams params, inout CloudSampleResult res);

void SampleCloudWithCheck(CloudSampleParams params, inout CloudSampleResult res) {
#if 1
	SampleCloud(params, res);
	res.sigma_t = clamp(params.height01, 0, 1) == params.height01 ? res.sigma_t : 0;
#else
	if (clamp(height01, 0, 1) != height01)
		return;
	SampleCloud(params, res);
#endif
}

#endif