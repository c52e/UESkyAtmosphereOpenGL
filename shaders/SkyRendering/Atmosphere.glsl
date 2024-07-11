#ifndef _ATMOSPHERE_GLSL
#define _ATMOSPHERE_GLSL

/*
This code contains portions of source code from
https://github.com/ebruneton/precomputed_atmospheric_scattering/blob/master/atmosphere/functions.glsl
Copyright (c) 2017 Eric Bruneton
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:
1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holders nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.
*/

#include "Common.glsl"

#define STORE_OPTICAL_LENGTH_IN_TRANSMITTANCE_TEXTURE 1
#define COMBINED_SCATTERING_TEXTURES 0

layout(std140, binding = 0) uniform AtmosphereBufferData{
    vec3 solar_illuminance;
    float sun_angular_radius;

    vec3 rayleigh_scattering;
    float inv_rayleigh_exponential_distribution;

    vec3 mie_scattering;
    float inv_mie_exponential_distribution;

    vec3 mie_absorption;
    float ozone_center_altitude;

    vec3 ozone_absorption;
    float inv_ozone_width;

    vec3 ground_albedo;
    float mie_phase_g;

    vec3 _atmosphere_padding;
    float multiscattering_mask;

    float bottom_radius;
    float top_radius;
    float transmittance_steps;
    float multiscattering_steps;
};

const float kMuSMin = cos(120.0 / 180.0 * PI);

const int SCATTERING_TEXTURE_R_SIZE = 32;
const int SCATTERING_TEXTURE_MU_SIZE = 128;
const int SCATTERING_TEXTURE_MU_S_SIZE = 32;
#define SCATTERING_TEXTURE_NU_SIZE 8

const int SCATTERING_TEXTURE_WIDTH =
    SCATTERING_TEXTURE_NU_SIZE * SCATTERING_TEXTURE_MU_S_SIZE;
const int SCATTERING_TEXTURE_HEIGHT = SCATTERING_TEXTURE_MU_SIZE;
const int SCATTERING_TEXTURE_DEPTH = SCATTERING_TEXTURE_R_SIZE;

#define TRANSMITTANCE_FORMAT rgba16f
#define MULTISCATTERING_FORMAT rgba16f
#define SCATTERING_FORMAT rgba16f
#define IRRADIANCE_FORMAT rgba16f

// https://ebruneton.github.io/precomputed_atmospheric_scattering/
// https://github.com/sebh/UnrealEngineSkyAtmosphere

float ClampCosine(float mu) {
    return clamp(mu, -1.0, 1.0);
}

float ClampDistance(float d) {
    return max(d, 0.0);
}

float ClampRadius(float r) {
    return clamp(r, bottom_radius, top_radius);
}

float SafeSqrt(float a) {
    return sqrt(max(a, 0.0));
}

float GetUnitRangeFromTextureCoord(float u, int texture_size) {
    return (u - 0.5 / float(texture_size)) / (1.0 - 1.0 / float(texture_size));
}

float GetTextureCoordFromUnitRange(float x, int texture_size) {
    return 0.5 / float(texture_size) + x * (1.0 - 1.0 / float(texture_size));
}

vec2 GetTextureCoordFromUnitRange(vec2 xy, ivec2 texture_size) {
    return 0.5 / vec2(texture_size) + xy * (1.0 - 1.0 / vec2(texture_size));
}

bool RayIntersectsGround(float r, float mu) {
    return mu < 0.0 && r * r * (mu * mu - 1.0) + bottom_radius * bottom_radius >= 0.0;
}

float DistanceToTopAtmosphereBoundary(float r, float mu) {
    float discriminant = r * r * (mu * mu - 1.0) + top_radius * top_radius;
    return ClampDistance(-r * mu + SafeSqrt(discriminant));
}

float DistanceToBottomAtmosphereBoundary(float r, float mu) {
    float discriminant = r * r * (mu * mu - 1.0) + bottom_radius * bottom_radius;
    return ClampDistance(-r * mu - SafeSqrt(discriminant));
}

float DistanceToNearestAtmosphereBoundary(
        float r, float mu, bool ray_r_mu_intersects_ground) {
    return ray_r_mu_intersects_ground ? 
        DistanceToBottomAtmosphereBoundary(r, mu) :
        DistanceToTopAtmosphereBoundary(r, mu);
}

void GetRMuFromTransmittanceTextureIndex(ivec2 index, ivec2 size, out float r, out float mu) {
    vec2 uv = vec2(index) / vec2(size - 1);
    float x_mu = uv.x;
    float x_r = uv.y;
    // Distance to top atmosphere boundary for a horizontal ray at ground level.
    float H = sqrt(top_radius * top_radius - bottom_radius * bottom_radius);
    // Distance to the horizon, from which we can compute r:
    float rho = H * x_r;
    r = sqrt(rho * rho + bottom_radius * bottom_radius);
    // Distance to the top atmosphere boundary for the ray (r,mu), and its minimum
    // and maximum values over all mu - obtained for (r,1) and (r,mu_horizon) -
    // from which we can recover mu:
    float d_min = top_radius - r;
    float d_max = rho + H;
    float d = d_min + x_mu * (d_max - d_min);
    mu = d == 0.0 ? 1.0 : (H * H - rho * rho - d * d) / (2.0 * r * d);
    mu = ClampCosine(mu);
}

vec2 GetTransmittanceTextureUvFromRMu(ivec2 size, float r, float mu) {
    // Distance to top atmosphere boundary for a horizontal ray at ground level.
    float H = sqrt(top_radius * top_radius - bottom_radius * bottom_radius);
    // Distance to the horizon.
    float rho = SafeSqrt(r * r - bottom_radius * bottom_radius);
    // Distance to the top atmosphere boundary for the ray (r,mu), and its minimum
    // and maximum values over all mu - obtained for (r,1) and (r,mu_horizon).
    float d = DistanceToTopAtmosphereBoundary(r, mu);
    float d_min = top_radius - r;
    float d_max = rho + H;
    float x_mu = (d - d_min) / (d_max - d_min);
    float x_r = rho / H;
    return GetTextureCoordFromUnitRange(vec2(x_mu, x_r), size);
}

vec3 GetTransmittanceToTopAtmosphereBoundary(sampler2D transmittance_texture, float r, float mu) {
    vec2 uv = GetTransmittanceTextureUvFromRMu(textureSize(transmittance_texture, 0), r, mu);
#if STORE_OPTICAL_LENGTH_IN_TRANSMITTANCE_TEXTURE
    return exp(-texture(transmittance_texture, uv).rgb);
#else
    return texture(transmittance_texture, uv).rgb;
#endif
}

#if STORE_OPTICAL_LENGTH_IN_TRANSMITTANCE_TEXTURE
vec3 GetOpticalLengthToTopAtmosphereBoundary(sampler2D transmittance_texture, float r, float mu) {
    vec2 uv = GetTransmittanceTextureUvFromRMu(textureSize(transmittance_texture, 0), r, mu);
    return texture(transmittance_texture, uv).rgb;
}
#endif

vec3 GetTransmittance(sampler2D transmittance_texture,
        float r, float mu, float d, bool ray_r_mu_intersects_ground) {

    float r_d = ClampRadius(sqrt(d * d + 2.0 * r * mu * d + r * r));
    float mu_d = ClampCosine((r * mu + d) / r_d);

#if STORE_OPTICAL_LENGTH_IN_TRANSMITTANCE_TEXTURE
    if (ray_r_mu_intersects_ground) {
        return min(
            exp(-GetOpticalLengthToTopAtmosphereBoundary(
                    transmittance_texture, r_d, -mu_d) +
                GetOpticalLengthToTopAtmosphereBoundary(
                    transmittance_texture, r, -mu)),
            vec3(1.0));
    }
    else {
        return min(
            exp(-GetOpticalLengthToTopAtmosphereBoundary(
                transmittance_texture, r, mu) +
            GetOpticalLengthToTopAtmosphereBoundary(
                transmittance_texture, r_d, mu_d)),
            vec3(1.0));
    }
#else
    if (ray_r_mu_intersects_ground) {
        return min(
            GetTransmittanceToTopAtmosphereBoundary(
                transmittance_texture, r_d, -mu_d) /
            GetTransmittanceToTopAtmosphereBoundary(
                transmittance_texture, r, -mu),
            vec3(1.0));
    }
    else {
        return min(
            GetTransmittanceToTopAtmosphereBoundary(
                transmittance_texture, r, mu) /
            GetTransmittanceToTopAtmosphereBoundary(
                transmittance_texture, r_d, mu_d),
            vec3(1.0));
    }
#endif
}

vec3 GetSunVisibility(sampler2D transmittance_texture, float r, float mu_s) {
    float sin_theta_h = bottom_radius / r;
    float cos_theta_h = -sqrt(max(1.0 - sin_theta_h * sin_theta_h, 0.0));
    return GetTransmittanceToTopAtmosphereBoundary(transmittance_texture, r, mu_s) *
        smoothstep(-sin_theta_h * sun_angular_radius,
            sin_theta_h * sun_angular_radius,
            mu_s - cos_theta_h);
}

vec3 GetExtinctionCoefficient(float altitude) {
    vec3 rayleigh_extinction = rayleigh_scattering *
        clamp(exp(-altitude * inv_rayleigh_exponential_distribution), 0, 1);

    vec3 mie_extinction = (mie_scattering + mie_absorption) *
        clamp(exp(-altitude * inv_mie_exponential_distribution), 0, 1);

    vec3 ozone_extinction = ozone_absorption *
        max(0, altitude < ozone_center_altitude ?
            1 + (altitude - ozone_center_altitude) * inv_ozone_width :
            1 - (altitude - ozone_center_altitude) * inv_ozone_width);

    return rayleigh_extinction + mie_extinction + ozone_extinction;
}

float IsotropicPhaseFunction() {
    return 1.0 / (4.0 * PI);
}

float RayleighPhaseFunction(float cos_theta) {
    float k = 3.0 / (16.0 * PI);
    return k * (1.0 + cos_theta * cos_theta);
}

float CornetteShanks(float cos_theta, float g) {
    float k = 3.0 / (8.0 * PI) * (1.0 - g * g) / (2.0 + g * g);
    return k * (1.0 + cos_theta * cos_theta) / pow(1.0 + g * g - 2.0 * g * cos_theta, 1.5);
}

float HenyeyGreenstein(float cos_theta, float g) {
    float a = 1.0 - g * g;
    float b = 1.0 + g * g - 2.0 * g * cos_theta;
    b *= sqrt(b);
    return (0.25 * INV_PI) * a / b;
}

float MiePhaseFunction(float cos_theta, float g) {
#if 0
    return CornetteShanks(cos_theta, g);
#else
    return HenyeyGreenstein(cos_theta, g);
#endif
}

void GetScatteringCoefficient(float altitude, out vec3 rayleigh, out vec3 mie) {
    rayleigh = rayleigh_scattering * clamp(exp(-altitude * inv_rayleigh_exponential_distribution), 0, 1);
    mie = mie_scattering * clamp(exp(-altitude * inv_mie_exponential_distribution), 0, 1);
}

void GetAltitudeMuSFromMultiscatteringTextureIndex(ivec2 index, ivec2 size, out float altitude, out float mu_s) {
    vec2 uv = vec2(index) / vec2(size - 1);
    float x_mu_s = uv.x;
    float x_altitude = uv.y;
    altitude = x_altitude * (top_radius - bottom_radius);
    mu_s = x_mu_s * 2.0 - 1.0;
}

vec2 GetMultiscatteringTextureUvFromRMuS(ivec2 size, float r, float mu_s) {
    float x_mu_s = mu_s * 0.5 + 0.5;
    float x_r = (r - bottom_radius) / (top_radius - bottom_radius);
    return GetTextureCoordFromUnitRange(vec2(x_mu_s, x_r), size);
}

vec3 GetMultiscatteringContribution(sampler2D multiscattering_texture, float r, float mu_s) {
    vec2 uv = GetMultiscatteringTextureUvFromRMuS(textureSize(multiscattering_texture, 0), r, mu_s);
    return texture(multiscattering_texture, uv).rgb;
}

float GetVisibilityFromShadowMap(sampler2DShadow shadow_map, mat4 light_view_projection, vec3 position) {
    vec3 xyz = MulDivW(light_view_projection, position);
    xyz.xy = xyz.xy * 0.5 + 0.5;
    vec2 shadow_map_coord = xyz.xy;
    float depth = xyz.z;
    if (depth >= 1.0) return 1.0;
    return texture(shadow_map, vec3(shadow_map_coord, depth));
}

float GetVisibilityFromMoonShadow(float sun_moon_angular_distance,
        float sun_angular_radius, float moon_angular_radius) {
    // 可以用一个预计算的二维纹理加速
    float max_radius = sun_angular_radius + moon_angular_radius;
    float min_radius = abs(sun_angular_radius - moon_angular_radius);
    float sun_r2 = sun_angular_radius * sun_angular_radius;
    float moon_r2 = moon_angular_radius * moon_angular_radius;
    if (sun_moon_angular_distance >= max_radius)
        return 1.0;
    if (sun_moon_angular_distance <= min_radius)
        return clamp((sun_r2 - moon_r2) / sun_r2, 0.0, 1.0);
    float distance2 = sun_moon_angular_distance * sun_moon_angular_distance;
    float cos_half_sun = (distance2 + sun_r2 - moon_r2) / (2 * sun_moon_angular_distance * sun_angular_radius);
    float cos_half_moon = (distance2 + moon_r2 - sun_r2) / (2 * sun_moon_angular_distance * moon_angular_radius);
    float half_sun = acos(cos_half_sun);
    float half_moon = acos(cos_half_moon);
    float triangle_h = sun_angular_radius * sqrt(1 - cos_half_sun * cos_half_sun);
    float area_total = (PI - half_sun) * sun_r2 + (PI - half_moon) * moon_r2 + triangle_h * sun_moon_angular_distance;
    float area_uncovered = area_total - PI * moon_r2;
    return area_uncovered / (PI * sun_r2);
}

float GetVisibilityFromMoonShadow(vec3 moon_vector, float moon_radius, vec3 sun_direction) {
    float inv_moon_distance = inversesqrt(dot(moon_vector, moon_vector));
    vec3 moon_direction = moon_vector * inv_moon_distance;
    float moon_angular_radius = asin(clamp(moon_radius * inv_moon_distance, -1, 1));
    return GetVisibilityFromMoonShadow(acos(clamp(dot(sun_direction, moon_direction), -1, 1)),
        sun_angular_radius, moon_angular_radius);
}

vec3 ComputeScatteredLuminance(sampler2D transmittance_texture
#ifndef MULTISCATTERING_COMPUTE_PROGRAM
    , sampler2D multiscattering_texture
#if VOLUMETRIC_LIGHT_ENABLE
    , sampler2DShadow shadow_map
    , mat4 light_view_projection
#endif
#if MOON_SHADOW_ENABLE
    , vec3 moon_position
    , float moon_radius
#endif
    , float start_i
#endif
    , vec3 earth_center, vec3 start_position, vec3 view_direction, vec3 sun_direction,
    float marching_distance, float steps, out vec3 transmittance
#ifdef MULTISCATTERING_COMPUTE_PROGRAM
    , out vec3 L_f
#endif
) {
    float r = length(start_position - earth_center);
    float rmu = dot(view_direction, start_position - earth_center);
    float mu_s = dot(view_direction, sun_direction);

    const float SAMPLE_COUNT = steps;
    float dx = marching_distance / SAMPLE_COUNT;

    transmittance = vec3(1.0);
    vec3 luminance = vec3(0.0);
#ifdef MULTISCATTERING_COMPUTE_PROGRAM
    L_f = vec3(0.0);
    float start_i = 0.5;
    float rayleigh_phase = IsotropicPhaseFunction();
    float mie_phase = IsotropicPhaseFunction();
#else
    float rayleigh_phase = RayleighPhaseFunction(mu_s);
    float mie_phase = MiePhaseFunction(mu_s, mie_phase_g);
#endif
    for (float i = start_i; i < SAMPLE_COUNT; ++i) {
        float d_i = i * dx;
        // Distance between the current sample point and the planet center.
        float r_i = sqrt(d_i * d_i + 2.0 * rmu * d_i + r * r);
        vec3 position_i = start_position + view_direction * d_i;
        float altitude_i = r_i - bottom_radius;

        vec3 rayleigh_scattering_i;
        vec3 mie_scattering_i;
        GetScatteringCoefficient(altitude_i, rayleigh_scattering_i, mie_scattering_i);
        vec3 scattering_i = rayleigh_scattering_i + mie_scattering_i;
        vec3 scattering_with_phase_i = rayleigh_scattering_i * rayleigh_phase + mie_scattering_i * mie_phase;

        vec3 extinction_i = GetExtinctionCoefficient(altitude_i);
        vec3 transmittance_i = exp(-extinction_i * dx);
        vec3 up_direction_i = normalize(position_i - earth_center);
        float mu_s_i = dot(sun_direction, up_direction_i);
        vec3 luminance_i = scattering_with_phase_i * GetSunVisibility(transmittance_texture, r_i, mu_s_i);
#ifndef MULTISCATTERING_COMPUTE_PROGRAM
#if VOLUMETRIC_LIGHT_ENABLE
        // Shadow Map是小范围阴影，不应影响多重散射的贡献
        luminance_i *= GetVisibilityFromShadowMap(shadow_map, light_view_projection, position_i);
#endif
        vec3 multiscattering_contribution = GetMultiscatteringContribution(
            multiscattering_texture, r_i, mu_s_i);
        luminance_i += multiscattering_mask * multiscattering_contribution * scattering_i;
#if MOON_SHADOW_ENABLE
        // 月亮阴影是大范围阴影，应该影响到多重散射的贡献
        luminance_i *= GetVisibilityFromMoonShadow(moon_position - position_i, moon_radius, sun_direction);
#endif
        luminance_i *= solar_illuminance;
#endif

        luminance += transmittance * (luminance_i - luminance_i * transmittance_i) / extinction_i;
#ifdef MULTISCATTERING_COMPUTE_PROGRAM
        L_f += transmittance * (scattering_i - scattering_i * transmittance_i) / extinction_i;
#endif
        transmittance *= transmittance_i;
    }
    return luminance;
}


vec4 GetScatteringTextureUvwzFromRMuMuSNu(
        float r, float mu, float mu_s, float nu, bool ray_r_mu_intersects_ground) {
    // Distance to top atmosphere boundary for a horizontal ray at ground level.
    float H = sqrt(top_radius * top_radius - bottom_radius * bottom_radius);
    // Distance to the horizon.
    float rho = SafeSqrt(r * r - bottom_radius * bottom_radius);
    float u_r = GetTextureCoordFromUnitRange(rho / H, SCATTERING_TEXTURE_R_SIZE);

    // Discriminant of the quadratic equation for the intersections of the ray
    // (r,mu) with the ground (see RayIntersectsGround).
    float r_mu = r * mu;
    float discriminant = r_mu * r_mu - r * r + bottom_radius * bottom_radius;
    float u_mu;
    if (ray_r_mu_intersects_ground) {
        // Distance to the ground for the ray (r,mu), and its minimum and maximum
        // values over all mu - obtained for (r,-1) and (r,mu_horizon).
        float d = -r_mu - SafeSqrt(discriminant);
        float d_min = r - bottom_radius;
        float d_max = rho;
        u_mu = 0.5 - 0.5 * GetTextureCoordFromUnitRange(d_max == d_min ? 0.0 :
            (d - d_min) / (d_max - d_min), SCATTERING_TEXTURE_MU_SIZE / 2);
    }
    else {
        // Distance to the top atmosphere boundary for the ray (r,mu), and its
        // minimum and maximum values over all mu - obtained for (r,1) and
        // (r,mu_horizon).
        float d = -r_mu + SafeSqrt(discriminant + H * H);
        float d_min = top_radius - r;
        float d_max = rho + H;
        u_mu = 0.5 + 0.5 * GetTextureCoordFromUnitRange(
            (d - d_min) / (d_max - d_min), SCATTERING_TEXTURE_MU_SIZE / 2);
    }

    float d = DistanceToTopAtmosphereBoundary(bottom_radius, mu_s);
    float d_min = top_radius - bottom_radius;
    float d_max = H;
    float a = (d - d_min) / (d_max - d_min);
    float D = DistanceToTopAtmosphereBoundary(bottom_radius, kMuSMin);
    float A = (D - d_min) / (d_max - d_min);
    // An ad-hoc function equal to 0 for mu_s = mu_s_min (because then d = D and
    // thus a = A), equal to 1 for mu_s = 1 (because then d = d_min and thus
    // a = 0), and with a large slope around mu_s = 0, to get more texture 
    // samples near the horizon.
    float u_mu_s = GetTextureCoordFromUnitRange(
        max(1.0 - a / A, 0.0) / (1.0 + a), SCATTERING_TEXTURE_MU_S_SIZE);

    float u_nu = (nu + 1.0) / 2.0;
    return vec4(u_nu, u_mu_s, u_mu, u_r);
}

vec3 GetExtrapolatedSingleMieScattering(vec4 scattering) {
    // Algebraically this can never be negative, but rounding errors can produce
    // that effect for sufficiently short view rays.
    if (scattering.r <= 0.0) {
        return vec3(0.0);
    }
    return scattering.rgb * scattering.a / scattering.r *
        (rayleigh_scattering.r / mie_scattering.r) *
        (mie_scattering / rayleigh_scattering);
}

vec3 GetCombinedScattering(sampler3D scattering_texture,
        sampler3D single_mie_scattering_texture,
        float r, float mu, float mu_s, float nu,
        bool ray_r_mu_intersects_ground,
        out vec3 single_mie_scattering) {
    vec4 uvwz = GetScatteringTextureUvwzFromRMuMuSNu(
        r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    float tex_coord_x = uvwz.x * float(SCATTERING_TEXTURE_NU_SIZE - 1);
    float tex_x = floor(tex_coord_x);
    float lerp = tex_coord_x - tex_x;
    vec3 uvw0 = vec3((tex_x + uvwz.y) / float(SCATTERING_TEXTURE_NU_SIZE),
        uvwz.z, uvwz.w);
    vec3 uvw1 = vec3((tex_x + 1.0 + uvwz.y) / float(SCATTERING_TEXTURE_NU_SIZE),
        uvwz.z, uvwz.w);
#if COMBINED_SCATTERING_TEXTURES
    vec4 combined_scattering =
#if SCATTERING_TEXTURE_NU_SIZE == 1
        texture(scattering_texture, uvw0);
#else
        texture(scattering_texture, uvw0) * (1.0 - lerp) +
        texture(scattering_texture, uvw1) * lerp;
#endif
    vec3 scattering = vec3(combined_scattering);
    single_mie_scattering = GetExtrapolatedSingleMieScattering(combined_scattering);
#else
    vec3 scattering = vec3(
#if SCATTERING_TEXTURE_NU_SIZE == 1
        texture(scattering_texture, uvw0));
#else
        texture(scattering_texture, uvw0) * (1.0 - lerp) +
        texture(scattering_texture, uvw1) * lerp);
#endif
    single_mie_scattering = vec3(
#if SCATTERING_TEXTURE_NU_SIZE == 1
        texture(single_mie_scattering_texture, uvw0));
#else
        texture(single_mie_scattering_texture, uvw0) * (1.0 - lerp) +
        texture(single_mie_scattering_texture, uvw1) * lerp);
#endif
#endif
    return scattering;
}

vec3 GetSkyRadiance(sampler2D transmittance_texture,
        sampler3D scattering_texture, sampler3D single_mie_scattering_texture,
        vec3 camera, vec3 view_ray, vec3 sun_direction, out vec3 transmittance) {
    // Compute the distance to the top atmosphere boundary along the view ray,
    // assuming the viewer is in space (or NaN if the view ray does not intersect
    // the atmosphere).
    float r = length(camera);
    float rmu = dot(camera, view_ray);

    // Compute the r, mu, mu_s and nu parameters needed for the texture lookups.
    float mu = rmu / r;
    float mu_s = dot(camera, sun_direction) / r;
    float nu = dot(view_ray, sun_direction);
    bool ray_r_mu_intersects_ground = RayIntersectsGround(r, mu);

    transmittance = ray_r_mu_intersects_ground ? vec3(0.0) :
        GetTransmittanceToTopAtmosphereBoundary(transmittance_texture, r, mu);
    vec3 single_mie_scattering;
    vec3 scattering;
    scattering = GetCombinedScattering(
        scattering_texture, single_mie_scattering_texture,
        r, mu, mu_s, nu, ray_r_mu_intersects_ground,
        single_mie_scattering);
    return scattering * RayleighPhaseFunction(nu) + single_mie_scattering *
        MiePhaseFunction(nu, mie_phase_g);
}

vec3 GetSkyRadianceToPoint(sampler2D transmittance_texture,
        sampler3D scattering_texture, sampler3D single_mie_scattering_texture,
        vec3 camera, vec3 view_ray, float d, vec3 sun_direction, out vec3 transmittance) {
    // Compute the distance to the top atmosphere boundary along the view ray,
    // assuming the viewer is in space (or NaN if the view ray does not intersect
    // the atmosphere).
    float r = length(camera);
    float rmu = dot(camera, view_ray);
    float distance_to_top_atmosphere_boundary = -rmu -
        sqrt(rmu * rmu - r * r + top_radius * top_radius);

    // Compute the r, mu, mu_s and nu parameters for the first texture lookup.
    float mu = rmu / r;
    float mu_s = dot(camera, sun_direction) / r;
    float nu = dot(view_ray, sun_direction);
    bool ray_r_mu_intersects_ground = RayIntersectsGround(r, mu);

    transmittance = GetTransmittance(transmittance_texture,
        r, mu, d, ray_r_mu_intersects_ground);

    vec3 single_mie_scattering;
    vec3 scattering = GetCombinedScattering(
        scattering_texture, single_mie_scattering_texture,
        r, mu, mu_s, nu, ray_r_mu_intersects_ground,
        single_mie_scattering);

    // Compute the r, mu, mu_s and nu parameters for the second texture lookup.
    float r_p = ClampRadius(sqrt(d * d + 2.0 * r * mu * d + r * r));
    float mu_p = (r * mu + d) / r_p;
    float mu_s_p = (r * mu_s + d * nu) / r_p;

    vec3 single_mie_scattering_p;
    vec3 scattering_p = GetCombinedScattering(
        scattering_texture, single_mie_scattering_texture,
        r_p, mu_p, mu_s_p, nu, ray_r_mu_intersects_ground,
        single_mie_scattering_p);

    // Combine the lookup results to get the scattering between camera and point.
    scattering = scattering - transmittance * scattering_p;
    single_mie_scattering =
        single_mie_scattering - transmittance * single_mie_scattering_p;
#if COMBINED_SCATTERING_TEXTURES
    single_mie_scattering = GetExtrapolatedSingleMieScattering(
        vec4(scattering, single_mie_scattering.r));
#endif

    // Hack to avoid rendering artifacts when the sun is below the horizon.
    single_mie_scattering = single_mie_scattering *
        smoothstep(float(0.0), float(0.01), mu_s);

    return scattering * RayleighPhaseFunction(nu) + single_mie_scattering *
        MiePhaseFunction(nu, mie_phase_g);
}


vec3 ComputeGroundLuminance(sampler2D transmittance_texture, vec3 earth_center, vec3 position, vec3 sun_direction) {
    vec3 up_direction = normalize(position - earth_center);
    float mu_s = dot(sun_direction, up_direction);
    vec3 solar_illuminance_at_ground = GetSunVisibility(transmittance_texture, bottom_radius, mu_s);
#ifndef MULTISCATTERING_COMPUTE_PROGRAM
    solar_illuminance_at_ground *= solar_illuminance;
#endif
    vec3 normal = normalize(position - earth_center);
    return INV_PI * clamp(dot(normal, sun_direction), 0, 1) * ground_albedo * solar_illuminance_at_ground;
}

void GetRMuMuSNuFromScatteringTextureUvwz(vec4 uvwz,
    out float r, out float mu, out float mu_s, out float nu, out bool ray_r_mu_intersects_ground) {
    // Distance to top atmosphere boundary for a horizontal ray at ground level.
    float H = sqrt(top_radius * top_radius - bottom_radius * bottom_radius);
    // Distance to the horizon.
    float rho =
        H * GetUnitRangeFromTextureCoord(uvwz.w, SCATTERING_TEXTURE_R_SIZE);
    r = sqrt(rho * rho + bottom_radius * bottom_radius);

    if (uvwz.z < 0.5) {
        // Distance to the ground for the ray (r,mu), and its minimum and maximum
        // values over all mu - obtained for (r,-1) and (r,mu_horizon) - from which
        // we can recover mu:
        float d_min = r - bottom_radius;
        float d_max = rho;
        float d = d_min + (d_max - d_min) * GetUnitRangeFromTextureCoord(
            1.0 - 2.0 * uvwz.z, SCATTERING_TEXTURE_MU_SIZE / 2);
        mu = d == 0.0 ? -1.0 : ClampCosine(-(rho * rho + d * d) / (2.0 * r * d));
        ray_r_mu_intersects_ground = true;
    }
    else {
        // Distance to the top atmosphere boundary for the ray (r,mu), and its
        // minimum and maximum values over all mu - obtained for (r,1) and
        // (r,mu_horizon) - from which we can recover mu:
        float d_min = top_radius - r;
        float d_max = rho + H;
        float d = d_min + (d_max - d_min) * GetUnitRangeFromTextureCoord(
            2.0 * uvwz.z - 1.0, SCATTERING_TEXTURE_MU_SIZE / 2);
        mu = d == 0.0 ? 1.0 : ClampCosine((H * H - rho * rho - d * d) / (2.0 * r * d));
        ray_r_mu_intersects_ground = false;
    }

    float x_mu_s =
        GetUnitRangeFromTextureCoord(uvwz.y, SCATTERING_TEXTURE_MU_S_SIZE);
    float d_min = top_radius - bottom_radius;
    float d_max = H;
    float D = DistanceToTopAtmosphereBoundary(bottom_radius, kMuSMin);
    float A = (D - d_min) / (d_max - d_min);
    float a = (A - x_mu_s * A) / (1.0 + x_mu_s * A);
    float d = d_min + min(a, A) * (d_max - d_min);
    mu_s = d == 0.0 ? 1.0 : ClampCosine((H * H - d * d) / (2.0 * bottom_radius * d));

    nu = ClampCosine(uvwz.x * 2.0 - 1.0);
}

void GetRMuMuSNuFromScatteringTextureFragCoord(
    vec3 frag_coord,
    out float r, out float mu, out float mu_s, out float nu,
    out bool ray_r_mu_intersects_ground) {
    const vec4 SCATTERING_TEXTURE_SIZE = vec4(
        SCATTERING_TEXTURE_NU_SIZE - 1,
        SCATTERING_TEXTURE_MU_S_SIZE,
        SCATTERING_TEXTURE_MU_SIZE,
        SCATTERING_TEXTURE_R_SIZE);
    float frag_coord_nu =
        floor(frag_coord.x / float(SCATTERING_TEXTURE_MU_S_SIZE));
    float frag_coord_mu_s =
        mod(frag_coord.x, float(SCATTERING_TEXTURE_MU_S_SIZE));
    vec4 uvwz =
        vec4(frag_coord_nu, frag_coord_mu_s, frag_coord.y, frag_coord.z) /
        SCATTERING_TEXTURE_SIZE;
    GetRMuMuSNuFromScatteringTextureUvwz(
        uvwz, r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    // Clamp nu to its valid range of values, given mu and mu_s.
    nu = clamp(nu, mu * mu_s - sqrt((1.0 - mu * mu) * (1.0 - mu_s * mu_s)),
        mu * mu_s + sqrt((1.0 - mu * mu) * (1.0 - mu_s * mu_s)));
}

void GetRMuSFromIrradianceTextureUv(vec2 frag_coord, ivec2 image_size, out float r, out float mu_s) {
    vec2 uv = frag_coord / image_size;
    float x_mu_s = GetUnitRangeFromTextureCoord(uv.x, image_size.x);
    float x_r = GetUnitRangeFromTextureCoord(uv.y, image_size.y);
    r = bottom_radius + x_r * (top_radius - bottom_radius);
    mu_s = ClampCosine(2.0 * x_mu_s - 1.0);
}

vec3 GetScattering(
    sampler3D scattering_texture,
    float r, float mu, float mu_s, float nu,
    bool ray_r_mu_intersects_ground) {
    vec4 uvwz = GetScatteringTextureUvwzFromRMuMuSNu(
        r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    float tex_coord_x = uvwz.x * float(SCATTERING_TEXTURE_NU_SIZE - 1);
    float tex_x = floor(tex_coord_x);
    float lerp = tex_coord_x - tex_x;
    vec3 uvw0 = vec3((tex_x + uvwz.y) / float(SCATTERING_TEXTURE_NU_SIZE),
        uvwz.z, uvwz.w);
    vec3 uvw1 = vec3((tex_x + 1.0 + uvwz.y) / float(SCATTERING_TEXTURE_NU_SIZE),
        uvwz.z, uvwz.w);
    return vec3(texture(scattering_texture, uvw0) * (1.0 - lerp) +
        texture(scattering_texture, uvw1) * lerp);
}

vec3 GetScattering(
    sampler3D single_rayleigh_scattering_texture,
    sampler3D single_mie_scattering_texture,
    sampler3D multiple_scattering_texture,
    float r, float mu, float mu_s, float nu,
    bool ray_r_mu_intersects_ground,
    int scattering_order) {
    if (scattering_order == 1) {
        vec3 rayleigh = GetScattering(
            single_rayleigh_scattering_texture, r, mu, mu_s, nu,
            ray_r_mu_intersects_ground);
        vec3 mie = GetScattering(
            single_mie_scattering_texture, r, mu, mu_s, nu,
            ray_r_mu_intersects_ground);
        return rayleigh * RayleighPhaseFunction(nu) +
            mie * MiePhaseFunction(nu, mie_phase_g);
    }
    else {
        return GetScattering(
            multiple_scattering_texture, r, mu, mu_s, nu,
            ray_r_mu_intersects_ground);
    }
}
#endif
