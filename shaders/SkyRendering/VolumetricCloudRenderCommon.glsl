#include "VolumetricCloudCommon.glsl"
#include "VolumetricCloudShadowInterface.glsl"

layout(std140, binding = 2) uniform VolumetricCloudBufferData {
	float uSunIlluminanceScale;
	float uMaxRaymarchDistance;
	float uBaseStepSize;
    float uMaxVisibleDistance;

	vec3 uEnvColorScale;
	float uShadowDistance;
    
    float uEnvMultiscatteringSigmaScale;
	float uEnvBottomVisibility;
	float uEnvSunHeightCurveExp;
    float pad1;
};

float SunCosTheta(vec3 pos) {
    vec3 up = normalize(vec3(pos.xy, pos.z + uEarthRadius));
    return clamp(dot(up, uSunDirection), 0, 1);
}

float SampleOpticalDepth(vec3 pos) {
    float optical_depth = 0.0;
    float inv_shadow_steps = 1.0 / SHADOW_STEPS;
    vec3 sample_vector = uShadowDistance * uSunDirection;
    float previous_t = 0.0;
    // UE4 Non-linear shadow sample distribution
    #pragma optionNV(unroll all)
    for (float t = inv_shadow_steps; t <= 1.0; t += inv_shadow_steps) {
        float current_t = t * t;
        float delta_t = current_t - previous_t;
        vec3 sample_pos = pos + sample_vector * (previous_t + 0.5 * delta_t);
        float sample_height01 = CalHeight01(sample_pos);
        
        CloudSampleResult res = DefaultCloudSampleResult();
        SampleCloudWithCheck(CloudSampleParams(sample_pos, sample_height01, uCameraPos, CLOUD_SAMPLE_SHADOW_RAY), res);
        optical_depth += res.sigma_t * uShadowDistance * delta_t;
        previous_t = current_t;
    }
    return optical_depth;
}
