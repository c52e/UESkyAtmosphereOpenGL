#include "Atmosphere.glsl"
#include "AtmosphereInterface.glsl"
#include "Shadow.glsl"
#include "VolumetricCloudShadowInterface.glsl"

layout(binding = 0) uniform sampler2D transmittance_texture;
layout(binding = 1) uniform sampler2D multiscattering_texture;
layout(binding = 2) uniform sampler2D depth_texture;
layout(binding = 3) uniform sampler2D albedo_texture;
layout(binding = 4) uniform sampler2D normal_texture;
layout(binding = 5) uniform sampler2D orm_texture;
layout(binding = 6) uniform sampler2DShadow shadow_map_texture;
layout(binding = 7) uniform sampler2D blue_noise;
layout(binding = 8) uniform sampler2D star_luminance;
layout(binding = 9) uniform sampler2D sky_view_luminance_texture;
layout(binding = 10) uniform sampler2D sky_view_transmittance_texture;
layout(binding = 11) uniform sampler3D aerial_perspective_luminance_texture;
layout(binding = 12) uniform sampler3D aerial_perspective_transmittance_texture;
layout(binding = 13) uniform sampler2D shadow_map_depth_sampler;
layout(binding = 14) uniform sampler3D cloud_shadow_froxel;
layout(binding = 15) uniform samplerCube prefiltered_radiance_texture;
layout(binding = 16) uniform sampler2D env_brdf_lut;
layout(binding = 17) uniform sampler2D pixel_visibility_tex;
layout(binding = 18) uniform sampler3D scattering_texture;
layout(binding = 19) uniform sampler3D single_mie_scattering_texture;

layout(std140, binding = 1) uniform AtmosphereRenderBufferData{
    vec3 sun_direction;
    float star_luminance_scale;
    vec3 earth_center;
    float camera_earth_center_distance;
    vec3 camera_position;
    float raymarching_steps;
    vec3 up_direction;
    float sky_view_lut_steps;
    vec3 right_direction;
    float aerial_perspective_lut_steps;
    vec3 front_direction;
    float aerial_perspective_lut_max_distance;
    vec3 moon_position;
    float moon_radius;

    mat4 inv_view_projection;
    mat4 light_view_projection;

    float padding00;
    float uInvShadowFroxelMaxDistance;
    float blocker_kernel_size_k;
    float pcss_size_k;
};

#define SKY_VIEW_FORMAT rgba16f
#define AERIAL_PERSPECTIVE_FORMAT rgba16f

layout(std140, binding = 2) uniform EnvRadianceSH {
    vec4 Llm[9];
};

vec3 ComputeRaymarchingStartPositionAndChangeDistance(vec3 view_direction, inout float marching_distance) {
    vec3 start_position = camera_position;
    bool in_space = camera_earth_center_distance > top_radius;
    if (in_space) {
        float r = camera_earth_center_distance;
        float mu = dot(view_direction, up_direction);
        float near_distance;
        if (FromSpaceIntersectTopAtmosphereBoundary(r, mu, near_distance)) {
            start_position += near_distance * view_direction;
            marching_distance -= near_distance;
        }
        else {
            marching_distance = 0;
        }
    }
    return start_position;
}

vec3 ComputeScatteredLuminance(
    bool is_sky, vec3 start_position, vec3 view_direction, float marching_distance, float steps, out vec3 transmittance
) {
#if USE_SCATTERING_LUT
    if (is_sky) {
        return GetSkyRadiance(transmittance_texture, scattering_texture, single_mie_scattering_texture,
            start_position - earth_center, view_direction, sun_direction, transmittance);
    }
    else {
        return GetSkyRadianceToPoint(transmittance_texture, scattering_texture, single_mie_scattering_texture,
            start_position - earth_center, view_direction, marching_distance, sun_direction, transmittance);
    }
#else
#if DITHER_SAMPLE_POINT_ENABLE
    float start_i = texelFetch(blue_noise, ivec2(gl_GlobalInvocationID.xy) & 0x3f, 0).x;
#else
    float start_i = 0.5;
#endif
    return ComputeScatteredLuminance(transmittance_texture, multiscattering_texture,
#if VOLUMETRIC_LIGHT_ENABLE
        shadow_map_texture, light_view_projection,
#endif
#if MOON_SHADOW_ENABLE
        moon_position, moon_radius,
#endif
        start_i, earth_center, start_position, view_direction,
        sun_direction, marching_distance, steps, transmittance);
#endif
}

float GetHorizonDownAngleFromR(float r) {
    float tangent_point_distance = sqrt(r * r - bottom_radius * bottom_radius);
    float cos_horizon_down = tangent_point_distance / r;
    return acos(cos_horizon_down);
}

void GetCosLatLonFromSkyViewTextureIndex(ivec2 index, float r, out float cos_lat, out float cos_lon) {
    vec2 uv = vec2(index) / vec2(SKY_VIEW_LUT_SIZE - 1);
    float x_cos_lon = uv.x;
    float x_cos_lat = uv.y;

    float horizon_down_angle = GetHorizonDownAngleFromR(r);
    float horizon_up_angle = PI - horizon_down_angle;

    float lat;
    if (x_cos_lat < 0.5f) {
        float coord = 1.0 - 2.0 * x_cos_lat;
        coord = 1.0 - coord * coord;
        lat = horizon_up_angle * coord;
    }
    else {
        float coord = x_cos_lat * 2.0 - 1.0;
        coord *= coord;
        lat = horizon_up_angle + horizon_down_angle * coord;
    }
    cos_lat = cos(lat);
    cos_lon = -(x_cos_lon * x_cos_lon * 2.0 - 1.0);
}

vec3 GetViewDirectionFromCosLatLon(float cos_lat, float cos_lon) {
    float sin_lat = clamp(sqrt(1 - cos_lat * cos_lat), 0, 1);
    float sin_lon = clamp(sqrt(1 - cos_lon * cos_lon), 0, 1);
    vec3 view_direction = up_direction * cos_lat +
        front_direction * (sin_lat * cos_lon) +
        right_direction * (sin_lat * sin_lon);
    return view_direction;
}

vec2 GetSkyViewTextureUvFromCosLatLon(float r, float cos_lat, float cos_lon) {
    float horizon_down_angle = GetHorizonDownAngleFromR(r);
    float horizon_up_angle = PI - horizon_down_angle;

    float lat = acos(cos_lat);
    float x_cos_lat;
    if (lat < horizon_up_angle) {
        float coord = lat / horizon_up_angle;
        coord = sqrt(1 - coord);
        x_cos_lat = 0.5 - 0.5 * coord;
    }
    else {
        float coord = (lat - horizon_up_angle) / horizon_down_angle;
        coord = sqrt(coord);
        x_cos_lat = coord * 0.5 + 0.5;
    }

    float x_cos_lon = sqrt(0.5 - 0.5 * cos_lon);
    return GetTextureCoordFromUnitRange(vec2(x_cos_lon, x_cos_lat), SKY_VIEW_LUT_SIZE);
}

void GetCosLatLonFromViewDirection(vec3 view_direction, out float cos_lat, out float cos_lon) {
    cos_lat = dot(up_direction, view_direction);
    vec3 lon_direction = view_direction - up_direction * cos_lat;
    float lon_direction_length2 = dot(lon_direction, lon_direction);
    if (lon_direction_length2 == 0) {
        cos_lon = 1;
    }
    else {
        lon_direction *= inversesqrt(lon_direction_length2);
        cos_lon = dot(lon_direction, front_direction);
    }
}

#ifdef SKY_VIEW_COMPUTE_PROGRAM
layout(binding = 0, SKY_VIEW_FORMAT) uniform image2D luminance_image;
layout(binding = 1, SKY_VIEW_FORMAT) uniform image2D transmittance_image;

void main() {
    float r = camera_earth_center_distance;
    float cos_lat, cos_lon;
    GetCosLatLonFromSkyViewTextureIndex(ivec2(gl_GlobalInvocationID.xy), r, cos_lat, cos_lon);
    vec3 view_direction = GetViewDirectionFromCosLatLon(cos_lat, cos_lon);

    float mu = cos_lat;
    bool intersect_bottom = RayIntersectsGround(r, mu);
    float marching_distance = DistanceToNearestAtmosphereBoundary(r, mu, intersect_bottom);

    vec3 start_position = ComputeRaymarchingStartPositionAndChangeDistance(view_direction, marching_distance);
    vec3 transmittance = vec3(1);
    vec3 luminance = vec3(0);
    if (marching_distance > 0) {
        luminance = ComputeScatteredLuminance(
            true, start_position, view_direction, marching_distance, sky_view_lut_steps, transmittance);
    }
    imageStore(luminance_image, ivec2(gl_GlobalInvocationID.xy), vec4(luminance, dot(transmittance, vec3(1.0/3.0))));
    imageStore(transmittance_image, ivec2(gl_GlobalInvocationID.xy), vec4(transmittance, 0.0));
}

#endif


float GetMarchingDistanceFromAerialPerspectiveTextureIndex(ivec3 index,
        vec3 camera_position, mat4 inv_view_projection, out vec3 view_direction) {
    vec3 uvw = vec3(index) / vec3(AERIAL_PERSPECTIVE_LUT_SIZE - 1);
    view_direction = UvToDir(inv_view_projection, uvw.xy);
    return uvw.z * uvw.z * aerial_perspective_lut_max_distance;
}

#ifdef AERIAL_PERSPECTIVE_COMPUTE_PROGRAM
layout(binding = 0, AERIAL_PERSPECTIVE_FORMAT) uniform image3D luminance_image;
layout(binding = 1, AERIAL_PERSPECTIVE_FORMAT) uniform image3D transmittance_image;

void main() {
    vec3 view_direction;
    float marching_distance = GetMarchingDistanceFromAerialPerspectiveTextureIndex(
        ivec3(gl_GlobalInvocationID.xyz), camera_position, inv_view_projection, view_direction);

    float r = camera_earth_center_distance;
    float mu = dot(view_direction, up_direction);
    bool intersect_bottom = RayIntersectsGround(r, mu);
    float max_marching_distance = intersect_bottom ?
#if 1
        marching_distance  // Better when camera is near ground
#else
        DistanceToBottomAtmosphereBoundary(r, mu) // Better when camera is near space
#endif
        : DistanceToTopAtmosphereBoundary(r, mu);
    vec3 start_position = ComputeRaymarchingStartPositionAndChangeDistance(view_direction, max_marching_distance);
    marching_distance = min(marching_distance, max_marching_distance);

    vec3 transmittance = vec3(1);
    vec3 luminance = vec3(0);
    if (marching_distance > 0) {
        luminance = ComputeScatteredLuminance(
            false, start_position, view_direction, marching_distance, aerial_perspective_lut_steps, transmittance);
    }
    imageStore(luminance_image, ivec3(gl_GlobalInvocationID.xyz), vec4(luminance, dot(transmittance, vec3(1.0 / 3.0))));
    imageStore(transmittance_image, ivec3(gl_GlobalInvocationID.xyz), vec4(transmittance, 0.0));
}

#endif


#ifdef ENVIRONMENT_LUMINANCE_COMPUTE_PROGRAM
layout(binding = 0, rgba16f) uniform imageCube env_luminance_image;

void main() {
    int index = int(gl_GlobalInvocationID.z);
    vec2 face_uv = (vec2(gl_GlobalInvocationID.xy) + vec2(0.5)) / vec2(imageSize(env_luminance_image));
    face_uv.y = 1.0 - face_uv.y;
    vec3 view_direction = ConvertCubUvToDir(index, face_uv);
    float cos_lat, cos_lon;
    GetCosLatLonFromViewDirection(view_direction, cos_lat, cos_lon);
    float r = camera_earth_center_distance;
    vec2 uv = GetSkyViewTextureUvFromCosLatLon(r, cos_lat, cos_lon);
    vec3 luminance = texture(sky_view_luminance_texture, uv).rgb;
    vec3 transmittance = texture(sky_view_transmittance_texture, uv).rgb;

    float mu = cos_lat;
    if (RayIntersectsGround(r, mu)) {
        float marching_distance = DistanceToBottomAtmosphereBoundary(r, mu);
        vec3 ground_position = camera_position + view_direction * marching_distance;
        luminance += ComputeGroundLuminance(transmittance_texture, earth_center, ground_position, sun_direction) * transmittance;
    }

    imageStore(env_luminance_image, ivec3(gl_GlobalInvocationID.xyz), vec4(luminance, 0.0));
}

#endif

#ifdef ATMOSPHERE_RENDER_FRAGMENT_SHADER

#include "../Base/BRDF.glsl"

layout(binding = 0, rgba16f) uniform image2D hdr_image;

float SampleVisibilityFromShadowMap(vec3 position) {
#if PCSS_ENABLE
    float visibility = PCSS(shadow_map_texture, shadow_map_depth_sampler, blocker_kernel_size_k, pcss_size_k, light_view_projection, position);
#else
    float visibility = GetVisibilityFromShadowMap(shadow_map_texture, light_view_projection, position);
#endif
    visibility = min(visibility, texelFetch(pixel_visibility_tex, ivec2(gl_GlobalInvocationID.xy), 0).r);
    return visibility;
}

vec3 ComputeObjectLuminance(vec3 position, vec3 view_direction) {
    float r = length(position - earth_center);
    vec3 object_up_direction = normalize(position - earth_center);
    float mu_s = dot(sun_direction, object_up_direction);

    vec3 sun_visibility;
    if (r > top_radius) {
        float near_distance;
        if (FromSpaceIntersectTopAtmosphereBoundary(r, mu_s, near_distance)) {
            position += near_distance * sun_direction;
            r = length(position - earth_center);
            object_up_direction = normalize(position - earth_center);
            mu_s = dot(sun_direction, object_up_direction);
            sun_visibility = GetSunVisibility(transmittance_texture, r, mu_s);
        }
        else {
            sun_visibility = vec3(1);
        }
    }
    else {
        sun_visibility = GetSunVisibility(transmittance_texture, r, mu_s);
    }

    vec3 solar_illuminance_at_object = solar_illuminance * sun_visibility;

    MaterialData material_data = LoadMeterialData(albedo_texture, normal_texture, orm_texture, ivec2(gl_GlobalInvocationID.xy));

    vec3 N = material_data.normal;
    vec3 L = sun_direction;
    vec3 V = -view_direction;
    vec3 H = normalize(L + V);
    float NdotL = clamp(dot(N, L), 0, 1);
    float NdotV = clamp(dot(N, V), 0, 1);
    float NdotH = clamp(dot(N, H), 0, 1);
    float HdotV = clamp(dot(H, V), 0, 1);
    vec3 brdf = BRDF(NdotL, NdotV, NdotH, HdotV, material_data);

    float shadow_visibility = 1.0;
    if (NdotL > 0) {
        shadow_visibility = SampleVisibilityFromShadowMap(position);
#if MOON_SHADOW_ENABLE
        shadow_visibility *= GetVisibilityFromMoonShadow(moon_position - position, moon_radius, sun_direction);
#endif
    }
    vec3 direct_lumiance = brdf * solar_illuminance_at_object * (NdotL * shadow_visibility);
    vec3 ambient_lumiance = GetAmbient(env_brdf_lut, ROUGHNESS_COUNT - 1, prefiltered_radiance_texture, V, material_data, Llm);
    float ambient_fade = clamp(10.0 - 0.1 * length(position - camera_position), 0, 1);
    return direct_lumiance + ambient_lumiance * ambient_fade;
}

vec3 GetStarLuminance(sampler2D star_luminance, vec3 view_direction) {
    float theta = acos(view_direction.y);
    float phi = atan(view_direction.x, view_direction.z);
    vec2 coord = vec2(INV_PI * 0.5 * phi + 0.5, 1.0 - theta * INV_PI);
    return star_luminance_scale * texture(star_luminance, coord).rgb;
}

void main() {
    ivec2 kImageSize = imageSize(hdr_image);
    if (gl_GlobalInvocationID.x >= kImageSize.x || gl_GlobalInvocationID.y >= kImageSize.y)
        return;
    vec2 vTexCoord = vec2(vec2(gl_GlobalInvocationID.xy) + 0.5) / vec2(kImageSize);

    float depth = texelFetch(depth_texture, ivec2(gl_GlobalInvocationID.xy), 0).x;
    bool intersect_object = depth != 0.0;
    vec3 fragment_position = UvDepthToPos(inv_view_projection, vTexCoord, depth);
    vec3 view_direction = UvToDir(inv_view_projection, vTexCoord);
    float r = camera_earth_center_distance;
    float mu = dot(view_direction, up_direction);

    float marching_distance = DistanceToTopAtmosphereBoundary(r, mu);
    if (intersect_object) {
        float object_distance = length(fragment_position - camera_position);
        marching_distance = min(marching_distance, object_distance);
    }

    vec3 start_position = ComputeRaymarchingStartPositionAndChangeDistance(view_direction, marching_distance);
    vec3 transmittance = vec3(1);
    vec3 luminance = vec3(0);
    if (marching_distance > 0) {
#if USE_SKY_VIEW_LUT
        if (!intersect_object) {
            float cos_lat, cos_lon;
            GetCosLatLonFromViewDirection(view_direction, cos_lat, cos_lon);
            vec2 uv = GetSkyViewTextureUvFromCosLatLon(r, cos_lat, cos_lon);
            luminance = texture(sky_view_luminance_texture, uv).rgb;
            transmittance = texture(sky_view_transmittance_texture, uv).rgb;
        }
        else
#endif
#if USE_AERIAL_PERSPECTIVE_LUT
        if (intersect_object) {
            vec3 uvw = GetAerialPerspectiveTextureUvwFromTexCoordDistance(
                vTexCoord, marching_distance, aerial_perspective_lut_max_distance, AERIAL_PERSPECTIVE_LUT_SIZE);
            luminance = texture(aerial_perspective_luminance_texture, uvw).rgb;
            transmittance = texture(aerial_perspective_transmittance_texture, uvw).rgb;
        }
        else
#endif
            luminance = ComputeScatteredLuminance(
                !intersect_object, start_position, view_direction, marching_distance, raymarching_steps, transmittance);

    }
    luminance *= SampleRayScatterVisibility(cloud_shadow_froxel, vTexCoord, marching_distance, uInvShadowFroxelMaxDistance);

    if (intersect_object) {
        luminance += transmittance * ComputeObjectLuminance(fragment_position, view_direction);
    }
    else if (dot(view_direction, sun_direction) >= cos(sun_angular_radius)) {
        // https://media.contentapi.ea.com/content/dam/eacom/frostbite/files/s2016-pbs-frostbite-sky-clouds-new.pdf
        // B Sun limb darkening astro-physical models
        vec3 u = vec3(1.0, 1.0, 1.0);
        vec3 a = vec3(0.397, 0.503, 0.652);

        float cos_view_sun = dot(view_direction, sun_direction);
        float sin_view_sun = sqrt(1.0 - cos_view_sun * cos_view_sun);
        
        float center_to_edge = clamp(sin_view_sun / sin(sun_angular_radius), 0.0, 1.0);
        float mu = sqrt(1.0 - center_to_edge * center_to_edge);

        vec3 factor = 1.0 - u * (1.0 - pow(vec3(mu), a));

        vec3 solar_illuminance_at_eye = solar_illuminance * transmittance;
        luminance += solar_illuminance_at_eye / (PI * sun_angular_radius * sun_angular_radius) * factor;
    } else {
        luminance += transmittance * GetStarLuminance(star_luminance, view_direction);
    }

    imageStore(hdr_image, ivec2(gl_GlobalInvocationID.xy), vec4(luminance, 1.0));
}

#endif
