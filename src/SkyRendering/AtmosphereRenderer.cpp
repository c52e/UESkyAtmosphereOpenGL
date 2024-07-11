#include "AtmosphereRenderer.h"

#include <array>
#include <sstream>

#include <glm/gtc/type_ptr.hpp>

#include "Utils.h"
#include "Textures.h"
#include "Samplers.h"
#include "VolumetricCloud.h"
#include "ImageLoader.h"

constexpr GLsizei kSkyViewTextureWidth = 128;
constexpr GLsizei kSkyViewTextureHeight = 128;
constexpr GLenum kSkyViewTextureInternalFormat = GL_RGBA16F;

constexpr GLsizei kAerialPerspectiveTextureWidth = 32;
constexpr GLsizei kAerialPerspectiveTextureHeight = 32;
constexpr GLenum kAerialPerspectiveTextureInternalFormat = GL_RGBA16F;

constexpr GLsizei kEnvironmentLuminanceTextureWidth = 128;

struct AtmosphereRenderBufferData {
    glm::vec3 sun_direction;
    float star_luminance_scale;
    glm::vec3 earth_center;
    float camera_earth_center_distance;
    glm::vec3 camera_position;
    float raymarching_steps;
    glm::vec3 up_direction;
    float sky_view_lut_steps;
    glm::vec3 right_direction;
    float aerial_perspective_lut_steps;
    glm::vec3 front_direction;
    float aerial_perspective_lut_max_distance;
    glm::vec3 moon_position;
    float moon_radius;

    glm::mat4 inv_view_projection;
    glm::mat4 light_view_projection;

    float padding00;
    float uInvShadowFroxelMaxDistance;
    float blocker_kernel_size_k;
    float pcss_size_k;
};

static void AssignBufferData(const AtmosphereRenderParameters& parameters, const Earth& earth, AtmosphereRenderBufferData& data) {
    FromThetaPhiToDirection(glm::radians(parameters.sun_direction_theta), 
         glm::radians(parameters.sun_direction_phi), glm::value_ptr(data.sun_direction));
    data.star_luminance_scale = parameters.star_luminance_scale;
    data.camera_position = parameters.camera_position;
    data.inv_view_projection = glm::inverse(parameters.view_projection);
    data.light_view_projection = parameters.light_view_projection;
    data.raymarching_steps = parameters.raymarching_steps;
    data.sky_view_lut_steps = parameters.sky_view_lut_steps;
    data.aerial_perspective_lut_steps = parameters.aerial_perspective_lut_steps;
    data.aerial_perspective_lut_max_distance = parameters.aerial_perspective_lut_max_distance;
    data.moon_position = glm::vec3(earth.moon_model()[3]);
    data.moon_radius = earth.moon_status.radius;
    data.blocker_kernel_size_k = parameters.blocker_kernel_size_k;
    data.pcss_size_k = parameters.pcss_size_k;

    data.earth_center = earth.center();
    data.camera_earth_center_distance = glm::distance(data.camera_position, data.earth_center);
    data.up_direction = glm::normalize(data.camera_position - data.earth_center);
    data.right_direction = glm::cross(data.sun_direction, data.up_direction);
    float right_direction_length = glm::length(data.right_direction);
    if (right_direction_length != 0.0) {
        data.right_direction /= right_direction_length;
    }
    else {
        glm::vec3 sun_direction_biased;
        FromThetaPhiToDirection( glm::radians(parameters.sun_direction_theta + 90.f),
            glm::radians(parameters.sun_direction_phi), glm::value_ptr(sun_direction_biased));
        data.right_direction = glm::normalize(glm::cross(sun_direction_biased, data.up_direction));
    }
    data.front_direction = glm::cross(data.up_direction, data.right_direction);
}

AtmosphereRenderer::AtmosphereRenderer(const AtmosphereRenderInitParameters& init_parameters)
    : use_sky_view_lut_(init_parameters.use_sky_view_lut), use_aerial_perspective_lut_(init_parameters.use_aerial_perspective_lut)
    , aerial_perspective_lut_depth_(init_parameters.aerial_perspective_lut_depth){
    atmosphere_render_buffer_.Create();
    glNamedBufferStorage(atmosphere_render_buffer_.id(), sizeof(AtmosphereRenderBufferData), NULL, GL_DYNAMIC_STORAGE_BIT);

    auto generate_shader_header = [init_parameters] (const std::string& header, bool dither_sample_point_enable) {
        std::stringstream ss;
        ss << "#version 460\n"
            << header
            << "#define SKY_VIEW_LUT_SIZE ivec2(" << kSkyViewTextureWidth << "," << kSkyViewTextureHeight << ")\n"
            << "#define AERIAL_PERSPECTIVE_LUT_SIZE ivec3(" << kAerialPerspectiveTextureWidth
            << "," << kAerialPerspectiveTextureHeight
            << "," << init_parameters.aerial_perspective_lut_depth << ")\n"
            << "#define PCSS_ENABLE " << (init_parameters.pcss_enable ? "1\n" : "0\n")
            << "#define VOLUMETRIC_LIGHT_ENABLE " << (init_parameters.volumetric_light_enable ? "1\n" : "0\n")
            << "#define MOON_SHADOW_ENABLE " << (init_parameters.moon_shadow_enable ? "1\n" : "0\n")
            << "#define DITHER_SAMPLE_POINT_ENABLE " << (dither_sample_point_enable ? "1\n" : "0\n")
            << "#define USE_SKY_VIEW_LUT " << (init_parameters.use_sky_view_lut ? "1\n" : "0\n")
            << "#define USE_AERIAL_PERSPECTIVE_LUT " << (init_parameters.use_aerial_perspective_lut ? "1\n" : "0\n")
            << "#define ROUGHNESS_COUNT " << IBL::kRoughnessCount << "\n"
            << "#define USE_SCATTERING_LUT " << (init_parameters.use_scattering_lut ? "1\n" : "0\n");
        return ss.str();
    };

    sky_view_program_ = {
        "../shaders/SkyRendering/AtmosphereRenderer.glsl",
        {{8, 4}, {8, 8}, {16, 4}, {16, 8}, {16, 16}, {32, 8}, {32, 16}},
        [generate_shader_header, dither = init_parameters.sky_view_lut_dither_sample_point_enable](const std::string& src) {
            return generate_shader_header("#define SKY_VIEW_COMPUTE_PROGRAM\n", dither) + src;
        },
        "SKY_VIEW"
    };

    sky_view_luminance_texture_.Create(GL_TEXTURE_2D);
    glTextureStorage2D(sky_view_luminance_texture_.id(), 1, kSkyViewTextureInternalFormat,
        kSkyViewTextureWidth, kSkyViewTextureHeight);

    sky_view_transmittance_texture_.Create(GL_TEXTURE_2D);
    glTextureStorage2D(sky_view_transmittance_texture_.id(), 1, kSkyViewTextureInternalFormat,
        kSkyViewTextureWidth, kSkyViewTextureHeight);

    aerial_perspective_program_ = {
        "../shaders/SkyRendering/AtmosphereRenderer.glsl",
        {{8, 4, 1}, {8, 8, 1}, {4, 4, 2}, {4, 4, 4}, {8, 4, 2}, {8, 4, 4}},
        [generate_shader_header, dither = init_parameters.aerial_perspective_lut_dither_sample_point_enable](const std::string& src) {
            return generate_shader_header("#define AERIAL_PERSPECTIVE_COMPUTE_PROGRAM\n", dither) + src;
        },
        "AERIAL_PERSPECTIVE"
    };

    aerial_perspective_luminance_texture_.Create(GL_TEXTURE_3D);
    glTextureStorage3D(aerial_perspective_luminance_texture_.id(), 1, kAerialPerspectiveTextureInternalFormat,
        kAerialPerspectiveTextureWidth, kAerialPerspectiveTextureHeight, aerial_perspective_lut_depth_);

    aerial_perspective_transmittance_texture_.Create(GL_TEXTURE_3D);
    glTextureStorage3D(aerial_perspective_transmittance_texture_.id(), 1, kAerialPerspectiveTextureInternalFormat,
        kAerialPerspectiveTextureWidth, kAerialPerspectiveTextureHeight, aerial_perspective_lut_depth_);

    environment_luminance_program_ = {
        "../shaders/SkyRendering/AtmosphereRenderer.glsl",
        {{8, 4, 1}, {8, 8, 1}, {16, 4, 1}, {16, 8, 1}},
        [generate_shader_header](const std::string& src) {
            return generate_shader_header("#define ENVIRONMENT_LUMINANCE_COMPUTE_PROGRAM\n", true) + src;
        },
    };

    glEnable(GL_TEXTURE_CUBE_MAP_SEAMLESS);
    environment_luminance_texture_.Create(GL_TEXTURE_CUBE_MAP);
    constexpr auto w = kEnvironmentLuminanceTextureWidth;
    glTextureStorage2D(environment_luminance_texture_.id(), GetMipmapLevels(w, w), GL_RGBA16F, w, w);

    render_program_ = {
        "../shaders/SkyRendering/AtmosphereRenderer.glsl",
        {{8, 8}, {8, 4}, {16, 4}, {16, 8}, {16, 16}, {32, 8}, {32, 16}},
        [generate_shader_header, dither = init_parameters.raymarching_dither_sample_point_enable](const std::string& src) {
            return generate_shader_header("#define ATMOSPHERE_RENDER_FRAGMENT_SHADER\n", dither) + src;
        },
        "VIEWPORT_RENDER"
    };

    froxel_accumulate_program_ = {
        "../shaders/SkyRendering/ShadowFroxelAccumulate.comp",
        {{8, 4, 1}, {8, 8, 1}, {16, 4, 1}, {16, 8, 1}},
    };
}

void AtmosphereRenderer::Render(const Earth& earth, const VolumetricCloud& volumetric_cloud, const AtmosphereRenderParameters& parameters) {
    AtmosphereRenderBufferData atmosphere_render_buffer_data_;
    AssignBufferData(parameters, earth, atmosphere_render_buffer_data_);

    auto cloud_shadow_froxel = volumetric_cloud.GetShadowFroxel();
    atmosphere_render_buffer_data_.uInvShadowFroxelMaxDistance = cloud_shadow_froxel.inv_max_distance;

    {
        PERF_MARKER("ShadowFroxelAccumulate")
        GLBindImageTextures({ cloud_shadow_froxel.shadow_froxel });
        glUseProgram(froxel_accumulate_program_.id());
        froxel_accumulate_program_.Dispatch(glm::vec2(cloud_shadow_froxel.dim));
    }

    sun_direction_ = atmosphere_render_buffer_data_.sun_direction;
    aerial_perspective_lut_max_distance_ = atmosphere_render_buffer_data_.aerial_perspective_lut_max_distance;

    glNamedBufferSubData(atmosphere_render_buffer_.id(), 0, sizeof(atmosphere_render_buffer_data_), &atmosphere_render_buffer_data_);
    glBindBufferBase(GL_UNIFORM_BUFFER, 1, atmosphere_render_buffer_.id());

    glBindBufferBase(GL_UNIFORM_BUFFER, 2, ibl_.env_radiance_sh_buffer());

    auto bind_textures = [this, &earth, &parameters, &cloud_shadow_froxel]() {
        GLBindTextures({ earth.atmosphere().transmittance_texture(),
                    earth.atmosphere().multiscattering_texture(),
                    parameters.gbuffer->depth(),
                    parameters.gbuffer->albedo(),
                    parameters.gbuffer->normal(),
                    parameters.gbuffer->orm(),
                    parameters.shadow_map_texture,
                    Textures::Instance().blue_noise(),
                    Textures::Instance().star_luminance(),
                    sky_view_luminance_texture_.id(),
                    sky_view_transmittance_texture_.id(),
                    aerial_perspective_luminance_texture_.id(),
                    aerial_perspective_transmittance_texture_.id(),
                    parameters.shadow_map_texture,
                    cloud_shadow_froxel.shadow_froxel,
                    ibl_.prefiltered_radiance(),
                    Textures::Instance().env_brdf_lut(),
                    parameters.gbuffer->pixel_visibility(),
                    earth.atmosphere().scattering_texture(),
                    earth.atmosphere().single_mie_scattering_texture(), });
        GLBindSamplers({ Samplers::GetLinearNoMipmapClampToEdge(),
                    Samplers::GetLinearNoMipmapClampToEdge(),
                    0u,
                    Samplers::GetLinearNoMipmapClampToEdge(),
                    Samplers::GetLinearNoMipmapClampToEdge(),
                    Samplers::GetLinearNoMipmapClampToEdge(),
                    Samplers::GetShadowMapSampler(),
                    0u,
                    Samplers::GetLinearNoMipmapClampToEdge(),
                    Samplers::GetLinearNoMipmapClampToEdge(),
                    Samplers::GetLinearNoMipmapClampToEdge(),
                    Samplers::GetLinearNoMipmapClampToEdge(),
                    Samplers::GetLinearNoMipmapClampToEdge(),
                    Samplers::GetNearestClampToEdge(),
                    cloud_shadow_froxel.sampler,
                    Samplers::GetAnisotropySampler(Samplers::Wrap::CLAMP_TO_EDGE),
                    Samplers::GetLinearNoMipmapClampToEdge(),
                    0u,
                    Samplers::GetLinearNoMipmapClampToEdge(),
                    Samplers::GetLinearNoMipmapClampToEdge() });
    };
    bind_textures();

    if (true) { // For environment luminance texture
        PERF_MARKER("SkyViewLut")
        GLBindImageTextures({ sky_view_luminance_texture_.id(), sky_view_transmittance_texture_.id() });
        glUseProgram(sky_view_program_.id());
        sky_view_program_.Dispatch({kSkyViewTextureWidth, kSkyViewTextureHeight});
    }

    if (true) { // For Volumetric Cloud
        PERF_MARKER("AerialPerspective")
        GLBindImageTextures({ aerial_perspective_luminance_texture_.id(), aerial_perspective_transmittance_texture_.id() });
        glUseProgram(aerial_perspective_program_.id());
        aerial_perspective_program_.Dispatch({ kAerialPerspectiveTextureWidth, kAerialPerspectiveTextureHeight, aerial_perspective_lut_depth_ });
    }
    glMemoryBarrier(GL_TEXTURE_FETCH_BARRIER_BIT);
    {
        PERF_MARKER("EnvironmentLuminance")
        GLBindImageTextures({ environment_luminance_texture_.id() });
        glUseProgram(environment_luminance_program_.id());
        environment_luminance_program_.Dispatch({ kEnvironmentLuminanceTextureWidth, kEnvironmentLuminanceTextureWidth, 6 });
    }
    glMemoryBarrier(GL_TEXTURE_FETCH_BARRIER_BIT);
    glGenerateTextureMipmap(environment_luminance_texture_.id());

    ibl_.Precompute(environment_luminance_texture_.id());
    bind_textures();
    {
        PERF_MARKER("Render")
        GLBindImageTextures({ parameters.hdrbuffer->hdr_texture() });
        glUseProgram(render_program_.id());
        render_program_.Dispatch(parameters.gbuffer->dim());
    }
    glMemoryBarrier(GL_TEXTURE_FETCH_BARRIER_BIT | GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
}
