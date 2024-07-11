#version 460
#include "Atmosphere.glsl"
#include "Terrain.glsl"

layout(std140, binding = 1) uniform EarthBufferData{
    mat4 view_projection;
    mat4 inv_view_projection;
    vec3 camera_position;
    float camera_earth_center_distance;
    vec3 earth_center;
    float padding;
    vec3 up_direction;
    float padding1;
    vec3 light_dir;
    float padding2;
};

layout(std140, binding = 2) uniform TerrainBufferData{
    float uTerrainMaxAltitude;
    float uMaxVisibleDistance;
    float uAltitudeLodK;
    float uNormalLodK;
    float uMarchingThreshholdTangent;
    float uNormalDeltaTangent;
    vec2 pad;
    vec3 uDirtAlbedo;
    float uMarchingStepScale;
    vec3 uGrassAlbedo;
    float uMaxSteps;
    vec3 uSnowAlbedo;
    float uTerrainFrequency;
};

float TerrainYLod(vec2 p, float lod) {
    return TerrainY(p, lod, uTerrainFrequency, uTerrainMaxAltitude, bottom_radius, earth_center);
}

float TerrainYAltitude(vec2 p, float eye_dist) {
    float lod = min(-log2(eye_dist * uAltitudeLodK), 11.0f);
    return TerrainYLod(p, lod);
}

float TerrainYNormal(vec2 p, float eye_dist) {
    float lod = min(-log2(eye_dist * uNormalLodK), 15.0f);
    return TerrainYLod(p, lod);
}

float CalcVisibility(vec3 start) {
    float visibility = 1.0;
    float k_soft = light_dir.y / tan(sun_angular_radius);
    float t = uTerrainMaxAltitude * 0.005;
    for (int i = 0; i < 100; i++) {
        vec3 p = start + t * light_dir;
        float ydiff = p.y - TerrainYAltitude(p.xz, length(p - camera_position));
        visibility = min(visibility, saturate(k_soft * ydiff / t));

        t += max(uMarchingStepScale * ydiff, 0.01 * uTerrainMaxAltitude);
        if (visibility < 1e-3 || length(p - earth_center) >= bottom_radius + uTerrainMaxAltitude)
            break;
    }
    return visibility;
}

#if SHADOW_PROGRAM
layout(binding = 0, r16) uniform image2D visibility_image;
layout(binding = 0) uniform sampler2D depth_texture;
layout(binding = 1) uniform sampler2D normal_texture;

void main() {
    ivec2 pos = ivec2(gl_GlobalInvocationID.xy);
    ivec2 image_size = imageSize(visibility_image);
    vec2 uv = (vec2(pos) + 0.5) / image_size;

    float depth = texelFetch(depth_texture, pos, 0).x;
    if (depth == 0.0)
        return;

    vec3 normal = texelFetch(normal_texture, pos, 0).xyz;
    if (dot(normal, light_dir) <= 0)
        return;

    vec3 frag_pos = UvDepthToPos(inv_view_projection, uv, depth);
    float dist = length(frag_pos - camera_position);


    float visibility = imageLoad(visibility_image, pos).r;
    visibility = min(visibility, CalcVisibility(frag_pos));
    imageStore(visibility_image, pos, vec4(visibility));
}

#elif SHADOW_FROXEL_PROGRAM

layout(binding = 0, r16) uniform image3D shadow_froxel_image;
layout(location = 0) uniform float uShadowFroxelMaxDistance;

void main() {
    ivec3 pos = ivec3(gl_GlobalInvocationID.xyz);
    ivec3 image_size = imageSize(shadow_froxel_image);
    vec2 uv = (vec2(pos.xy) + 0.5) / image_size.xy;
    vec3 dir = UvToDir(inv_view_projection, uv);

    float visibility = imageLoad(shadow_froxel_image, pos).r;
    float t = (0.5 + float(pos.z)) / image_size.z * uShadowFroxelMaxDistance;
    vec3 p = camera_position + t * dir;
    visibility = min(visibility, CalcVisibility(p));
    imageStore(shadow_froxel_image, pos, vec4(visibility));
}

#else

layout(binding = 0) uniform sampler2D depth_texture;
layout(binding = 1) uniform sampler2D earth_albedo;

in vec2 vTexCoord;
layout(location = 0) out vec4 Albedo;
layout(location = 1) out vec4 Normal;
layout(location = 2) out vec4 ORM;

vec3 GetEarthAlbedo(sampler2D earth_albedo, vec3 ground_position) {
    vec3 direction = normalize(ground_position - earth_center);
    float theta = acos(direction.y);
    float phi = atan(direction.x, direction.z);
    vec2 coord = vec2(INV_PI * 0.5 * phi + 0.5, 1.0 - theta * INV_PI);
#if 0
    vec3 color = texture(earth_albedo, coord).rgb;
#else
    // Make the Earth seamless
    vec2 dudxy1 = vec2(dFdx(coord.x), dFdy(coord.x));
    vec2 dudxy2 = vec2(dFdx(fract(coord.x + 0.5)), dFdy(fract(coord.x + 0.5)));
    vec2 dudxy = length(dudxy1) < length(dudxy2) ? dudxy1 : dudxy2;
    vec2 dvdxy = vec2(dFdx(coord.y), dFdy(coord.y));
    vec3 color = textureGrad(earth_albedo, coord, vec2(dudxy.x, dvdxy.x), vec2(dudxy.y, dvdxy.y)).rgb;
#endif
    return color;
}

float GetAltitude(vec3 p) {
    return length(p - earth_center) - bottom_radius;
}

// https://iquilezles.org/articles/terrainmarching/
float Raymarch(vec3 dir, float tmin, float tmax) {
    float t = tmin;
	for (float i = 0; i < uMaxSteps; i++) {
        vec3 p = camera_position + t * dir;
		float ydiff = p.y - TerrainYAltitude(p.xz, t);
		if (abs(ydiff) < uMarchingThreshholdTangent * t || t > tmax)
            break;
		t += uMarchingStepScale * ydiff;
	}
	return t;
}

vec3 GetNormal(vec3 p, float t) {
    float delta = uNormalDeltaTangent * t;
    return normalize(vec3(TerrainYNormal(vec2(p.x - delta, p.z), t) - TerrainYNormal(vec2(p.x + delta, p.z), t),
                        2.0 * delta,
                    TerrainYNormal(vec2(p.x, p.z - delta), t) - TerrainYNormal(vec2(p.x, p.z + delta), t)));
}

float IntersectTerrain(vec3 dir, out vec3 albedo, out vec3 normal) {
    float tmin = 0.001;
    float tmax = uMaxVisibleDistance;

    float R0 = bottom_radius - uTerrainMaxAltitude * 1.1;
    float R1 = bottom_radius + uTerrainMaxAltitude * 1.1;
    float r = camera_earth_center_distance;
    float mu = dot(normalize(camera_position - earth_center), dir);
    float discriminant0 = r * r * (mu * mu - 1.0) + R0 * R0;
    float discriminant1 = r * r * (mu * mu - 1.0) + R1 * R1;
    float sqrt_discriminant0 = sqrt(max(discriminant0, 0));
    float sqrt_discriminant1 = sqrt(max(discriminant1, 0));
    if (discriminant1 <= 0)
        return -1.0;
    if (discriminant0 >= 0.0 && mu < 0.0)
        tmax = min(tmax, -r * mu - sqrt_discriminant0);
    else
        tmax = min(tmax, -r * mu + sqrt_discriminant1);
    if (r >= R1)
        tmin = max(tmin, -r * mu - sqrt_discriminant1);

    float t = Raymarch(dir, tmin, tmax);
    if (t >= tmax)
        return -1.0;
			
	vec3 p = camera_position + t * dir;
    normal = GetNormal(p, t);
    
    float altitude = GetAltitude(p);
    //dirt
    albedo = uDirtAlbedo;
    {
        // grass
        float min_altitude = 0.0;
        float max_altitude = 0.66 * uTerrainMaxAltitude;
        float noise = 0.625 * (PerlinFBM(33 * uTerrainFrequency * p.xz, 4, 1) * 0.8 + 0.8);
        float altitude_blend = smoothstep(min_altitude, max_altitude, altitude + noise);
        float normal_blend = smoothstep(0.5 * altitude_blend, 0.9 * altitude_blend, normal.y);
        float blend = (1.0 - altitude_blend * 0.5) * normal_blend;
        albedo = mix(albedo, uGrassAlbedo, blend);
    }
	{
        // snow
        float min_altitude = 0.66 * uTerrainMaxAltitude;
        float max_altitude = 0.86 * uTerrainMaxAltitude;
        float noise = 0.625 * (PerlinFBM(3.333 * uTerrainFrequency * p.xz, 4, 1) * 0.8 + 0.8) ;
	    float altitude_blend = smoothstep(min_altitude, max_altitude, altitude + noise);
        float normal_blend = smoothstep(1.0 - 0.5 * altitude_blend, 1.0 - 0.1 * altitude_blend, normal.y);
        float dir_blend = smoothstep(-1.3 * normal_blend, -0.0 * normal_blend, normal.x);
        float blend = altitude_blend * normal_blend * dir_blend;
        albedo = mix(albedo, uSnowAlbedo, blend);
    }
	return t;
}


void main() {
    float depth = texelFetch(depth_texture, ivec2(gl_FragCoord.xy), 0).x;
    vec3 fragment_position = UvDepthToPos(inv_view_projection, vTexCoord, depth);
    float frag_dist = depth != 0.0 ? distance(fragment_position, camera_position) : 1e9;
    vec3 view_direction = UvToDir(inv_view_projection, vTexCoord);

    float metallic = 0.0;
    float roughness = 1.0;
    ORM = vec4(1.0, roughness, metallic, 1.0);

#if ENABLE_RAYMARCHED_TERRAIN
    float t = IntersectTerrain(view_direction, Albedo.rgb, Normal.xyz);
    if (t < 0 || t > frag_dist)
        discard;
    vec3 position = camera_position + view_direction * t;
    gl_FragDepth = MulDivW(view_projection, position).z;
#else
    float r = camera_earth_center_distance;
    float mu = dot(view_direction, up_direction);
    if (!RayIntersectsGround(r, mu))
        discard;
    float dist = DistanceToBottomAtmosphereBoundary(r, mu);
    if (dist >= frag_dist)
        discard;
    vec3 ground_position = camera_position + view_direction * dist;
    gl_FragDepth = MulDivW(view_projection, ground_position).z;
    vec3 albedo = GetEarthAlbedo(earth_albedo, ground_position);
    Albedo = vec4(albedo, 1.0);
    vec3 normal = normalize(ground_position - earth_center);
    Normal = vec4(normal, 1.0);
#endif
}

#endif
