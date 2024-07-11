#include "Atmosphere.h"

#include <array>

#include "Textures.h"
#include "Samplers.h"

constexpr glm::ivec2 kTransmittanceTextureDim{ 256, 64 };
constexpr GLenum kTransmittanceTextureInternalFormat = GL_RGBA16F;

constexpr glm::ivec2 kMultiscatteringTextureDim{ 32, 32 };
constexpr GLenum kMultiscatteringTextureInternalFormat = GL_RGBA16F;

constexpr int SCATTERING_TEXTURE_R_SIZE = 32;
constexpr int SCATTERING_TEXTURE_MU_SIZE = 128;
constexpr int SCATTERING_TEXTURE_MU_S_SIZE = 32;
constexpr int SCATTERING_TEXTURE_NU_SIZE = 8;
constexpr glm::ivec3 kScatteringTextureDim{
    SCATTERING_TEXTURE_NU_SIZE * SCATTERING_TEXTURE_MU_S_SIZE,
    SCATTERING_TEXTURE_MU_SIZE,
    SCATTERING_TEXTURE_R_SIZE
};
constexpr GLenum kScatteringTextureInternalFormat = GL_RGBA16F;

constexpr glm::ivec2 kIrradianceTextureDim{ 64, 16 };
constexpr GLenum kIrradianceTextureInternalFormat = GL_RGBA16F;

struct AtmosphereBufferData {
    glm::vec3 solar_illuminance;
    float sun_angular_radius;

    glm::vec3 rayleigh_scattering;
    float inv_rayleigh_exponential_distribution;

    glm::vec3 mie_scattering;
    float inv_mie_exponential_distribution;

    glm::vec3 mie_absorption;
    float ozone_center_altitude;

    glm::vec3 ozone_absorption;
    float inv_ozone_width;

    glm::vec3 ground_albedo;
    float mie_phase_g;

    glm::vec3 _atmosphere_padding;
    float multiscattering_mask;

    float bottom_radius;
    float top_radius;
    float transmittance_steps;
    float multiscattering_steps;
};

static void AssignBufferData(const AtmosphereParameters& parameters, AtmosphereBufferData& data) {
    data.solar_illuminance = parameters.solar_illuminance;
    data.sun_angular_radius = glm::radians(parameters.sun_angular_radius);

    data.rayleigh_scattering = parameters.rayleigh_scattering * parameters.rayleigh_scattering_scale;
    data.mie_scattering = parameters.mie_scattering * parameters.mie_scattering_scale;
    data.mie_absorption = parameters.mie_absorption * parameters.mie_absorption_scale;
    data.ozone_absorption = parameters.ozone_absorption * parameters.ozone_absorption_scale;

    data.inv_rayleigh_exponential_distribution = 1.0f / parameters.rayleigh_exponential_distribution;
    data.inv_mie_exponential_distribution = 1.0f / parameters.mie_exponential_distribution;
    data.ozone_center_altitude = parameters.ozone_center_altitude;
    data.inv_ozone_width = 1.0f / parameters.ozone_width;
    data.ground_albedo = parameters.ground_albedo;
    data.bottom_radius = parameters.bottom_radius;
    data.top_radius = parameters.bottom_radius + parameters.thickness;
    data.mie_phase_g = parameters.mie_phase_g;

    data.transmittance_steps = parameters.transmittance_steps;
    data.multiscattering_steps = parameters.multiscattering_steps;
    data.multiscattering_mask = parameters.multiscattering_mask;
}

Atmosphere::Atmosphere() {
    atmosphere_parameters_buffer_.Create();
    glNamedBufferStorage(atmosphere_parameters_buffer_.id(), sizeof(AtmosphereBufferData), NULL, GL_DYNAMIC_STORAGE_BIT);

    transmittance_texture_.Create(GL_TEXTURE_2D);
    glTextureStorage2D(transmittance_texture_.id(), 1, kTransmittanceTextureInternalFormat,
        kTransmittanceTextureDim.x, kTransmittanceTextureDim.y);

    transmittance_program_ = {
        "../shaders/SkyRendering/AtmosphereTransmittanceLut.comp",
        {{16, 8}, {8, 8}, {8, 4}}
    };

    multiscattering_texture_.Create(GL_TEXTURE_2D);
    glTextureStorage2D(multiscattering_texture_.id(), 1, kMultiscatteringTextureInternalFormat,
        kMultiscatteringTextureDim.x, kMultiscatteringTextureDim.y);

    multiscattering_program_ = []() {
        auto multiscattering_program_src = 
            ReadWithPreprocessor("../shaders/SkyRendering/AtmosphereMultiscatteringLut.comp");
        return GLProgram(multiscattering_program_src.c_str()); };

    scattering_texture_.Create(GL_TEXTURE_3D);
    glTextureStorage3D(scattering_texture_.id(), 1, kScatteringTextureInternalFormat,
        kScatteringTextureDim.x, kScatteringTextureDim.y, kScatteringTextureDim.z);
    single_mie_scattering_texture_.Create(GL_TEXTURE_3D);
    glTextureStorage3D(single_mie_scattering_texture_.id(), 1, kScatteringTextureInternalFormat,
        kScatteringTextureDim.x, kScatteringTextureDim.y, kScatteringTextureDim.z);

    scattering_program_ = {
        "../shaders/SkyRendering/AtmosphereScatteringLut.comp",
        {{8, 8}, {16, 8}, {8, 4}},
        [](const std::string& src) {return Replace(src, "USE_MULTISCATTERING_LUT", "1"); }
    };
    compute_single_scattering_ = {
        "../shaders/SkyRendering/AtmosphereScatteringLut.comp",
        {{8, 8}, {16, 8}, {8, 4}}
    };

    compute_direct_irradiance_ = {
        "../shaders/SkyRendering/AtmosphereDirectIrradiance.comp",
        {{8, 8}, {16, 8}, {8, 4}}
    };
    compute_scattering_density_ = {
        "../shaders/SkyRendering/AtmosphereScatteringDensity.comp",
        {{8, 8}, {16, 8}, {8, 4}}
    };
    compute_indirect_irradiance_ = {
        "../shaders/SkyRendering/AtmosphereIndirectIrradiance.comp",
        {{8, 8}, {16, 8}, {8, 4}}
    };
    compute_multiple_scattering_ = {
        "../shaders/SkyRendering/AtmosphereMultipleScattering.comp",
        {{8, 8}, {16, 8}, {8, 4}}
    };
    delta_irradiance_texture_.Create(GL_TEXTURE_2D);
    glTextureStorage2D(delta_irradiance_texture_.id(), 1, kIrradianceTextureInternalFormat,
        kIrradianceTextureDim.x, kIrradianceTextureDim.y);

    delta_multiple_scattering_texture_.Create(GL_TEXTURE_3D);
    glTextureStorage3D(delta_multiple_scattering_texture_.id(), 1, kScatteringTextureInternalFormat,
        kScatteringTextureDim.x, kScatteringTextureDim.y, kScatteringTextureDim.z);

    delta_scattering_density_texture_.Create(GL_TEXTURE_3D);
    glTextureStorage3D(delta_scattering_density_texture_.id(), 1, kScatteringTextureInternalFormat,
        kScatteringTextureDim.x, kScatteringTextureDim.y, kScatteringTextureDim.z);
}

void Atmosphere::UpdateLuts(const AtmosphereParameters& parameters) {
    PERF_MARKER("UpdateLuts")
    AtmosphereBufferData atmosphere_buffer_data_;
    AssignBufferData(parameters, atmosphere_buffer_data_);
    bool use_multiscattering_lut = parameters.num_scattering_orders == 0;

    glNamedBufferSubData(atmosphere_parameters_buffer_.id(), 0, sizeof(atmosphere_buffer_data_), &atmosphere_buffer_data_);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, atmosphere_parameters_buffer_.id());

    GLBindImageTextures({ transmittance_texture_.id() });
    {
        PERF_MARKER("Transmittance")
        glUseProgram(transmittance_program_.id());
        transmittance_program_.Dispatch(kTransmittanceTextureDim);
    }
    {
        PERF_MARKER("Multiscattering")
        GLBindImageTextures({ multiscattering_texture_.id() });
        GLBindTextures({ transmittance_texture() });
        GLBindSamplers({ Samplers::GetLinearNoMipmapClampToEdge() });
        glUseProgram(multiscattering_program_.id());
        glDispatchCompute(kMultiscatteringTextureDim.x, kMultiscatteringTextureDim.y, 1);
    }
    glMemoryBarrier(GL_TEXTURE_FETCH_BARRIER_BIT);
    const auto& compute_scattering = use_multiscattering_lut ? scattering_program_ : compute_single_scattering_;
    {
        PERF_MARKER("Scattering")
        GLBindImageTextures({ scattering_texture_.id(),
                            single_mie_scattering_texture_.id() });
        GLBindTextures({ transmittance_texture(),
                    multiscattering_texture() });
        GLBindSamplers({ Samplers::GetLinearNoMipmapClampToEdge(),
                    Samplers::GetLinearNoMipmapClampToEdge() });
        glUseProgram(compute_scattering.id());
        compute_scattering.Dispatch(kScatteringTextureDim);
    }
    glMemoryBarrier(GL_TEXTURE_FETCH_BARRIER_BIT);

    if (use_multiscattering_lut)
        return;
    {
        PERF_MARKER("DirectIrradiance");
        GLBindImageTextures({ delta_irradiance_texture_.id() });
        //GLBindTextures({ transmittance_texture() });
        //GLBindSamplers({ Samplers::GetLinearNoMipmapClampToEdge() });
        glUseProgram(compute_direct_irradiance_.id());
        compute_direct_irradiance_.Dispatch(kIrradianceTextureDim);
    }
    glMemoryBarrier(GL_TEXTURE_FETCH_BARRIER_BIT);
    for (int scattering_order = 2; scattering_order <= parameters.num_scattering_orders; ++scattering_order) {
        {
            PERF_MARKER("ScatteringDensity");
            GLBindImageTextures({ delta_scattering_density_texture_.id() });
            GLBindTextures({ transmittance_texture(),
                        scattering_texture(),
                        single_mie_scattering_texture(), 
                        delta_multiple_scattering_texture_.id(),
                        delta_irradiance_texture_.id() });
            GLBindSamplers({ Samplers::GetLinearNoMipmapClampToEdge(),
                        Samplers::GetLinearNoMipmapClampToEdge(),
                        Samplers::GetLinearNoMipmapClampToEdge(),
                        Samplers::GetLinearNoMipmapClampToEdge(),
                        Samplers::GetLinearNoMipmapClampToEdge() });
            glUseProgram(compute_scattering_density_.id());
            glUniform1i(0, scattering_order);
            compute_scattering_density_.Dispatch(kScatteringTextureDim);
        }
        glMemoryBarrier(GL_TEXTURE_FETCH_BARRIER_BIT);
        {
            PERF_MARKER("IndirectIrradiance");
            GLBindImageTextures({ delta_irradiance_texture_.id() });
            //GLBindTextures({});
            //GLBindSamplers({});
            glUseProgram(compute_indirect_irradiance_.id());
            glUniform1i(0, scattering_order - 1);
            compute_indirect_irradiance_.Dispatch(kIrradianceTextureDim);
        }
        glMemoryBarrier(GL_TEXTURE_FETCH_BARRIER_BIT);
        {
            PERF_MARKER("MultipleScattering");
            GLBindImageTextures({ delta_multiple_scattering_texture_.id(),
                            scattering_texture() });
            GLBindTextures({ delta_scattering_density_texture_.id() });
            GLBindSamplers({ Samplers::GetLinearNoMipmapClampToEdge() });
            glUseProgram(compute_multiple_scattering_.id());
            compute_multiple_scattering_.Dispatch(kScatteringTextureDim);
        }
        glMemoryBarrier(GL_TEXTURE_FETCH_BARRIER_BIT);
    }
}
