#ifndef _COMMON_GLSL
#define _COMMON_GLSL

#define PI 3.1415926535897932384626433832795
#define INV_PI (1.0 / PI)

#define saturate(x) clamp((x), 0, 1)

vec3 MulDivW(mat4 m, vec3 v) {
	vec4 xyzw = m * vec4(v, 1.0);
	return xyzw.xyz / xyzw.w;
}

vec3 UvDepthToPos(mat4 m, vec2 uv, float depth) {
    return MulDivW(m, vec3(uv * 2.0 - 1.0, depth));
}

// https://nlguillemot.wordpress.com/2016/12/07/reversed-z-in-opengl/
// reversed and infinity z
vec3 UvToDir(mat4 m, vec2 uv) {
    return normalize(vec3(m * vec4(uv * 2.0 - 1.0, 0.0, 1.0)));
}

// https://en.wikipedia.org/wiki/Cube_mapping
vec3 ConvertCubUvToDir(int index, vec2 uv) {
    // convert range 0 to 1 to -1 to 1
    float u = uv.x;
    float v = uv.y;
    float uc = 2.0 * u - 1.0;
    float vc = 2.0 * v - 1.0;
    vec3 dir = vec3(0);
    switch (index)
    {
    case 0: dir = vec3(1.0, vc, -uc); break;	// POSITIVE X
    case 1: dir = vec3(-1.0, vc, uc); break;	// NEGATIVE X
    case 2: dir = vec3(uc, 1.0, -vc); break;	// POSITIVE Y
    case 3: dir = vec3(uc, -1.0, vc); break;	// NEGATIVE Y
    case 4: dir = vec3(uc, vc, 1.0); break;	// POSITIVE Z
    case 5: dir = vec3(-uc, vc, -1.0); break;	// NEGATIVE Z
    }
    return normalize(dir);
}

void CreateOrthonormalBasis(vec3 N, out vec3 t0, out vec3 t1) {
#if 0
    const vec3 up = vec3(0.0, 1.0, 0.0);
    const vec3 right = vec3(1.0, 0.0, 0.0);
    if (abs(dot(N, up)) > 0.01) {
        t0 = normalize(cross(N, right));
        t1 = normalize(cross(N, t0));
    }
    else {
        t0 = normalize(cross(N, up));
        t1 = normalize(cross(N, t0));
    }
#else
    // https://graphics.pixar.com/library/OrthonormalB/paper.pdf
    float s = (N.z >= 0.0 ? 1.0 : -1.0);
    float a = -1.0 / (s + N.z);
    float b = N.x * N.y * a;
    t0 = vec3(1.0 + s * N.x * N.x * a, s * b, -s * N.x);
    t1 = vec3(b, s + N.y * N.y * a, -N.y);
#endif
}

// The Solid Angle of a Plane Triangle
// https://ieeexplore.ieee.org/document/4121581
float GetSolidAngle(vec3 R1, vec3 R2, vec3 R3) {
    float N = dot(R1, cross(R2, R3));
    float lenR1 = length(R1);
    float lenR2 = length(R2);
    float lenR3 = length(R3);
    float D = lenR1 * lenR2 * lenR3 + dot(R1, R2) * lenR3 + dot(R1, R3) * lenR2 + dot(R2, R3) * lenR1;
    return 2.0 * abs(atan(N, D));
}

vec3 Linear2sRGB(vec3 color) {
    return pow(color, vec3(1.0 / 2.2));
}

vec3 sRGB2Linear(vec3 color) {
    return pow(color, vec3(2.2));
}

#endif