
#include "../Base/Noise.glsl"

float PerlinFBM(vec2 p, float lod, uint seed) {
    float sum = 0.0;
    float weight = 0.5;
    while (lod > 0.0) {
        sum += saturate(lod--) * weight * (PerlinNoiseInfinity(p, seed));
		weight *= 0.5;
        p = mat2(0.8, -0.6, 0.6, 0.8) * (p * 2.0f);
    }
	return sum;
}

float TerrainY(vec2 p, float lod, float frequency, float max_altitude, float bottom_radius, vec3 earth_center) {
    float noise01 = PerlinFBM(frequency * p, lod, 0) + 0.5f;
    float altitude = max_altitude * noise01;
    altitude *= smoothstep(3.0f, 7.0f, length(p));
    vec2 v = p - vec2(earth_center.x, earth_center.z);
    return altitude + sqrt(max(bottom_radius * bottom_radius - dot(v, v), 0.0f)) + earth_center.y;
}
