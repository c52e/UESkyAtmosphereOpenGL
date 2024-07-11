#pragma once

#include <glm/glm.hpp>

#include "gl.hpp"
#include "GLReloadableProgram.h"
#include "Serialization.h"

struct AtmosphereParameters {
    // 角度单位：degree
    // 长度单位：km
    glm::vec3 solar_illuminance{ 1.f, 1.f, 1.f };
    float sun_angular_radius = 0.5334f * 0.5f;

    float bottom_radius = 6360.0f;
    float thickness = 6420.0f - 6360.0f;
    glm::vec3 ground_albedo{ 0.4f, 0.4f, 0.4f };

    float rayleigh_exponential_distribution = 8.0f;
    float rayleigh_scattering_scale = 0.0331f;
    glm::vec3 rayleigh_scattering{ 0.175287f, 0.409607f, 1.0f };

    float mie_exponential_distribution = 1.2f;
    float mie_phase_g = 0.8f;
    float mie_scattering_scale = 0.003996f;
    glm::vec3 mie_scattering{ 1.0f, 1.0f, 1.0f };
    float mie_absorption_scale = 0.000444f;
    glm::vec3 mie_absorption{ 1.0f, 1.0f, 1.0f };

    float ozone_center_altitude = 25.0f;
    float ozone_width = 15.0f;
    float ozone_absorption_scale = 0.001881f;
    glm::vec3 ozone_absorption{ 0.345561f, 1.0f , 0.045189f };

    float transmittance_steps = 40.0f;
    float multiscattering_steps = 30.0f;
    float multiscattering_mask = 1.0f; // 是否考虑多重散射[0,1]

    int num_scattering_orders = 0; // 0表示使用multiscattering lut
};

STRUCT_FIELD_DECLARATION_BEGIN(AtmosphereParameters, ISerializable)
    STRUCT_FIELD_DECLARE(solar_illuminance)
    STRUCT_FIELD_DECLARE(sun_angular_radius)

    STRUCT_FIELD_DECLARE(bottom_radius)
    STRUCT_FIELD_DECLARE(thickness)
    STRUCT_FIELD_DECLARE(ground_albedo)

    STRUCT_FIELD_DECLARE(rayleigh_exponential_distribution)
    STRUCT_FIELD_DECLARE(rayleigh_scattering_scale)
    STRUCT_FIELD_DECLARE(rayleigh_scattering)

    STRUCT_FIELD_DECLARE(mie_exponential_distribution)
    STRUCT_FIELD_DECLARE(mie_phase_g)
    STRUCT_FIELD_DECLARE(mie_scattering_scale)
    STRUCT_FIELD_DECLARE(mie_scattering)
    STRUCT_FIELD_DECLARE(mie_absorption_scale)
    STRUCT_FIELD_DECLARE(mie_absorption)

    STRUCT_FIELD_DECLARE(ozone_center_altitude)
    STRUCT_FIELD_DECLARE(ozone_width)
    STRUCT_FIELD_DECLARE(ozone_absorption_scale)
    STRUCT_FIELD_DECLARE(ozone_absorption)

    STRUCT_FIELD_DECLARE(transmittance_steps)
    STRUCT_FIELD_DECLARE(multiscattering_steps)
    STRUCT_FIELD_DECLARE(multiscattering_mask)

    STRUCT_FIELD_DECLARE(num_scattering_orders)
STRUCT_FIELD_DECLARATION_END()

class Atmosphere {
public:
    Atmosphere();

    void UpdateLuts(const AtmosphereParameters& parameters);

    GLuint transmittance_texture() const {
        return transmittance_texture_.id();
    }

    GLuint multiscattering_texture() const {
        return multiscattering_texture_.id();
    }

    GLuint scattering_texture() const {
        return scattering_texture_.id();
    }

    GLuint single_mie_scattering_texture() const {
        return single_mie_scattering_texture_.id();
    }

private:
    GLBuffer atmosphere_parameters_buffer_;

    GLTexture transmittance_texture_;
    GLReloadableComputeProgram transmittance_program_;

    GLTexture multiscattering_texture_;
    GLReloadableProgram multiscattering_program_;

    GLTexture scattering_texture_;
    GLTexture single_mie_scattering_texture_;
    GLReloadableComputeProgram scattering_program_;
    GLReloadableComputeProgram compute_single_scattering_;

    GLReloadableComputeProgram compute_direct_irradiance_;
    GLReloadableComputeProgram compute_scattering_density_;
    GLReloadableComputeProgram compute_indirect_irradiance_;
    GLReloadableComputeProgram compute_multiple_scattering_;

    GLTexture delta_irradiance_texture_;
    GLTexture delta_scattering_density_texture_;
    GLTexture delta_multiple_scattering_texture_;
};
