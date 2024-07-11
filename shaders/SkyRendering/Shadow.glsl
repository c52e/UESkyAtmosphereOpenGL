// games202 hw1
#ifndef _SHADOW_GLSL
#define _SHADOW_GLSL

#include "Common.glsl"

// Shadow map related variables
#define NUM_SAMPLES 25
#define NUM_RINGS 3

#define PI2 (PI * 2.0)

float rand_1to1(float x) {
    // -1 -1
    return fract(sin(x) * 10000.0);
}

float rand_2to1(vec2 uv) {
    // 0 - 1
    const float a = 12.9898, b = 78.233, c = 43758.5453;
    float dt = dot(uv.xy, vec2(a, b)), sn = mod(dt, PI);
    return fract(sin(sn) * c);
}

vec2 poissonDisk[NUM_SAMPLES];

void poissonDiskSamples(const in vec2 randomSeed) {

    float ANGLE_STEP = PI2 * float(NUM_RINGS) / float(NUM_SAMPLES);
    float INV_NUM_SAMPLES = 1.0 / float(NUM_SAMPLES);

    float angle = rand_2to1(randomSeed) * PI2;
    float radius = INV_NUM_SAMPLES;
    float radiusStep = radius;

    for (int i = 0; i < NUM_SAMPLES; i++) {
        poissonDisk[i] = vec2(cos(angle), sin(angle)) * pow(radius, 0.75);
        radius += radiusStep;
        angle += ANGLE_STEP;
    }
}

void uniformDiskSamples(const in vec2 randomSeed) {
    float randNum = rand_2to1(randomSeed);
    float sampleX = rand_1to1(randNum);
    float sampleY = rand_1to1(sampleX);

    float angle = sampleX * PI2;
    float radius = sqrt(sampleY);

    for (int i = 0; i < NUM_SAMPLES; i++) {
        poissonDisk[i] = vec2(radius * cos(angle), radius * sin(angle));

        sampleX = rand_1to1(sampleY);
        sampleY = rand_1to1(sampleX);

        angle = sampleX * PI2;
        radius = sqrt(sampleY);
    }
}

float FindBlocker(sampler2D shadow_map_depth_sampler, vec3 coords, float size) {
    float sum = 0.0;
    float cnt = 0.0;
    for (int i = 0; i < NUM_SAMPLES; ++i) {
        vec2 samplePos = poissonDisk[i] * size + coords.xy;
        float shadow_depth = texture2D(shadow_map_depth_sampler, samplePos).x;
        if (coords.z - shadow_depth > 0.0) {
            cnt += 1.0;
            sum += shadow_depth;
        }
    }
    return sum / max(cnt, 1e-5);
}

float Filtering(sampler2DShadow shadowMap, vec3 coords, float size) {
    float sum = 0.0;
    for (int i = 0; i < NUM_SAMPLES; ++i) {
        vec2 samplePos = poissonDisk[i] * size + coords.xy;
        sum += texture(shadowMap, vec3(samplePos, coords.z));
    }
    return sum / float(NUM_SAMPLES);
}

float PCSS(sampler2DShadow shadowMap, sampler2D shadow_map_depth_sampler, float blocker_kernel_size_k, float pcss_size_k
        , mat4 light_view_projection, vec3 position) {
    vec3 coords = MulDivW(light_view_projection, position);
    coords.xy = coords.xy * 0.5 + 0.5;
    if (coords.z >= 1.0) return 1.0;
    
    poissonDiskSamples(coords.xy);
    float kernelSizeApproximate = blocker_kernel_size_k * coords.z;
    float avgblockerDepth = FindBlocker(shadow_map_depth_sampler, coords, kernelSizeApproximate);
    float distanceToFragment = coords.z - avgblockerDepth;
    float penumbraSize = pcss_size_k * distanceToFragment;

    return Filtering(shadowMap, coords, penumbraSize);
}

#endif
