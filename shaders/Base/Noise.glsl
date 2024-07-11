#ifndef _NOISE_GLSL
#define _NOISE_GLSL

uint WangHash(uint seed) {
    seed = (seed ^ 61) ^ (seed >> 16);
    seed *= 9;
    seed = seed ^ (seed >> 4);
    seed *= 0x27d4eb2d;
    seed = seed ^ (seed >> 15);
    return seed;
}

// https://www.reedbeta.com/blog/hash-functions-for-gpu-rendering/
uint PCGHash(uint seed) {
    uint state = seed * 747796405u + 2891336453u;
    uint word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
    return (word >> 22u) ^ word;
}

const vec2 kPerlinGradients2D[4] = {
    {1, 1}, {-1, 1}, {1, -1}, {-1, -1}
};

vec2 GetPerlinGradients(uint i, uint j, uint seed) {
    return kPerlinGradients2D[WangHash(seed + WangHash(i + WangHash(j))) & 0x3];
}

// Perlin Noise 2d
float PerlinNoiseInfinity(vec2 p, uint seed) {
    uvec2 ij0 = uvec2(ivec2(floor(p))); // Casting a negative float to uint is undefined
    uvec2 ij1 = uvec2(ivec2(ceil(p)));
    uint i0 = ij0.x;
    uint j0 = ij0.y;
    uint i1 = ij1.x;
    uint j1 = ij1.y;

    vec2 t = p - floor(p);
    vec2 uv = t * t * t * (t * (t * 6.0f - 15.0f) + 10.0f);
    float u = uv.x;
    float v = uv.y;

    float x0 = t.x;
    float y0 = t.y;
    float x1 = t.x - 1.0f;
    float y1 = t.y - 1.0f;

    return mix(mix(dot(GetPerlinGradients(i0, j0, seed), vec2(x0, y0)),
                   dot(GetPerlinGradients(i1, j0, seed), vec2(x1, y0)), u),
               mix(dot(GetPerlinGradients(i0, j1, seed), vec2(x0, y1)),
                   dot(GetPerlinGradients(i1, j1, seed), vec2(x1, y1)), u), v);
}

const vec3 kPerlinGradients[16] = {
    {1, 1, 0}, {-1, 1, 0}, {1, -1, 0}, {-1, -1, 0},
    {1, 0, 1}, {-1, 0, 1}, {1, 0, -1}, {-1, 0, -1},
    {0, 1, 1}, {0, -1, 1}, {0, 1, -1}, {0, -1, -1},
    {1, 1, 0}, {-1, 1, 0}, {0, -1, 1}, {0, -1, -1}
};

vec3 GetPerlinGradients(uint i, uint j, uint k, uint seed) {
    return kPerlinGradients[WangHash(seed + WangHash(i + WangHash(j + WangHash(k)))) & 0xf];
}

// https://www.semanticscholar.org/paper/Improving-noise-Perlin/a6fd5071b73f542c79bd08d409c5f73de38dac5d
// Tiled Perlin Noise 3d
float PerlinNoise(vec3 p, uint freq, uint seed) {
    p *= float(freq);
    uvec3 ijk0 = uvec3(ivec3(floor(p))) % freq; // Casting a negative float to uint is undefined
    uvec3 ijk1 = uvec3(ivec3(ceil(p))) % freq;
    uint i0 = ijk0.x;
    uint j0 = ijk0.y;
    uint k0 = ijk0.z;
    uint i1 = ijk1.x;
    uint j1 = ijk1.y;
    uint k1 = ijk1.z;

    vec3 t = p - floor(p);
    vec3 uvw = t * t * t * (t * (t * 6.0f - 15.0f) + 10.0f);
    float u = uvw.x;
    float v = uvw.y;
    float w = uvw.z;

    float x0 = t.x;
    float y0 = t.y;
    float z0 = t.z;
    float x1 = t.x - 1.0f;
    float y1 = t.y - 1.0f;
    float z1 = t.z - 1.0f;

    return mix(mix(mix(dot(GetPerlinGradients(i0, j0, k0, seed), vec3(x0, y0, z0)),
                       dot(GetPerlinGradients(i1, j0, k0, seed), vec3(x1, y0, z0)), u),
                   mix(dot(GetPerlinGradients(i0, j1, k0, seed), vec3(x0, y1, z0)),
                       dot(GetPerlinGradients(i1, j1, k0, seed), vec3(x1, y1, z0)), u),v),
               mix(mix(dot(GetPerlinGradients(i0, j0, k1, seed), vec3(x0, y0, z1)),
                       dot(GetPerlinGradients(i1, j0, k1, seed), vec3(x1, y0, z1)), u),
                   mix(dot(GetPerlinGradients(i0, j1, k1, seed), vec3(x0, y1, z1)),
                       dot(GetPerlinGradients(i1, j1, k1, seed), vec3(x1, y1, z1)), u), v), w);
}

// Tiled Worley Noise 2d
float WorleyNoise(vec2 p, uint freq, uint seed) {
    p *= float(freq);
    uvec2 ij = uvec2(floor(p));
    p += float(freq);
    float min_dist = 1e10;
    for (uint di = freq - 1; di <= freq + 1; ++di) {
        for (uint dj = freq - 1; dj <= freq + 1; ++dj) {
            uvec2 grid = ij + uvec2(di, dj);
            uvec2 gridseed = grid % freq;
            uint rnd0 = WangHash(seed + WangHash(gridseed.x + WangHash(gridseed.y)));
            uint rnd1 = WangHash(seed + WangHash(gridseed.x + WangHash(gridseed.y + 1)));
            vec2 g = vec2(grid) + vec2(rnd0, rnd1) / 4294967296.0f;
            float dist = distance(p, g);
            min_dist = min(min_dist, dist);
        }
    }
    return min_dist;
}

// Tiled Worley Noise 3d
float WorleyNoise(vec3 p, uint freq, uint seed) {
    p *= float(freq);
    uvec3 ijk = uvec3(floor(p));
    p += float(freq);
    float min_dist = 1e10;
    for (uint di = freq - 1; di <= freq + 1; ++di) {
        for (uint dj = freq - 1; dj <= freq + 1; ++dj) {
            for (uint dk = freq - 1; dk <= freq + 1; ++dk) {
                uvec3 grid = ijk + uvec3(di, dj, dk);
                uvec3 gridseed = grid % freq;
                uint rnd0 = WangHash(seed + WangHash(gridseed.x + WangHash(gridseed.y + WangHash(gridseed.z))));
                uint rnd1 = WangHash(seed + WangHash(gridseed.x + WangHash(gridseed.y + WangHash(gridseed.z + 1))));
                uint rnd2 = WangHash(seed + WangHash(gridseed.x + WangHash(gridseed.y + WangHash(gridseed.z + 2))));
                vec3 g = vec3(grid) + vec3(rnd0, rnd1, rnd2) / 4294967296.0f;
                float dist = distance(p, g);
                min_dist = min(min_dist, dist);
            }
        }
    }
    return min_dist;
}

// UE4
uint ReverseBits32(uint bits) {
    bits = (bits << 16) | (bits >> 16);
    bits = ((bits & 0x00ff00ff) << 8) | ((bits & 0xff00ff00) >> 8);
    bits = ((bits & 0x0f0f0f0f) << 4) | ((bits & 0xf0f0f0f0) >> 4);
    bits = ((bits & 0x33333333) << 2) | ((bits & 0xcccccccc) >> 2);
    bits = ((bits & 0x55555555) << 1) | ((bits & 0xaaaaaaaa) >> 1);
    return bits;
}

vec2 Hammersley(uint Index, uint NumSamples, uvec2 Random) {
    float E1 = fract(float(Index) / NumSamples + float(Random.x & 0xffff) / (1 << 16));
    float E2 = float(ReverseBits32(Index) ^ Random.y) * 2.3283064365386963e-10f;
    return vec2(E1, E2);
}

#endif
