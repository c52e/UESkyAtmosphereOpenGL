#pragma once

#include <optional>

#include <glm/glm.hpp>

#include "gl.hpp"
#include "GLReloadableProgram.h"
#include "Camera.h"
#include "Atmosphere.h"
#include "GBuffer.h"

class Earth : public ISerializable {
public:
    AtmosphereParameters parameters;
    bool always_update_atmosphere = false;

    bool enable_raymarched_terrain = true;
    bool enable_terrain_shadow = true;

    struct TerrainParamaters : public ISerializable {
        float max_altitude = 3.0f;
        float slope = 0.39f;
        float max_visible_distance = 500.0f;
        float normal_delta_tangent = 0.001f;
        float altitude_fineness = 0.60485f;
        float normal_fineness = 0.70485f;
        float marching_fineness = 0.45596f;
        float physical_lod = 4.0f;
        glm::vec3 dirt_albedo{ 0.15f, 0.15f, 0.1f };
        glm::vec3 grass_albedo{ 0.05f, 0.15f, 0.03f };
        glm::vec3 snow_albedo{ 0.62f, 0.65f, 0.7f };

        FIELD_DECLARATION_BEGIN(ISerializable)
            FIELD_DECLARE(max_altitude)
            FIELD_DECLARE(slope)
            FIELD_DECLARE(max_visible_distance)
            FIELD_DECLARE(normal_delta_tangent)
            FIELD_DECLARE(altitude_fineness)
            FIELD_DECLARE(normal_fineness)
            FIELD_DECLARE(marching_fineness)
            FIELD_DECLARE(physical_lod)
            FIELD_DECLARE(dirt_albedo)
            FIELD_DECLARE(grass_albedo)
            FIELD_DECLARE(snow_albedo)
        FIELD_DECLARATION_END()
    } terrain;

    struct MoonStatus : public ISerializable {
        float direction_theta = 70.0f;
        float direction_phi = 150.0f;
        float distance = 384401.f / 10.f;
        float radius = 1737.f;

        FIELD_DECLARATION_BEGIN(ISerializable)
            FIELD_DECLARE(direction_theta)
            FIELD_DECLARE(direction_phi)
            FIELD_DECLARE(distance)
            FIELD_DECLARE(radius)
        FIELD_DECLARATION_END()
    } moon_status;

    FIELD_DECLARATION_BEGIN(ISerializable)
        FIELD_DECLARE(parameters)
        FIELD_DECLARE(moon_status)
        FIELD_DECLARE(terrain)
        FIELD_DECLARE(enable_raymarched_terrain)
        FIELD_DECLARE(enable_terrain_shadow)
        FIELD_DECLARE(always_update_atmosphere)
    FIELD_DECLARATION_END()

    Earth();

    void Update(const Camera& camera, glm::vec2 viewport, const glm::vec3& light_dir);

    void DrawTerrainGui();

    glm::vec3 CheckCollision(const glm::vec3& target_pos, float zNear);

    glm::vec3 center() const { return { 0.0f, -parameters.bottom_radius, 0.0f }; }

    void RenderToGBuffer(GLuint depth_texture);

    struct FroxelInfo {
        GLuint froxel;
        glm::ivec3 dim;
        float max_distance;
    };
    void RenderShadowBuffer(const GBuffer& gbuffer, const FroxelInfo& froxel_info);

    const Atmosphere& atmosphere() const { return atmosphere_; }

    glm::mat4 moon_model() const;

private:
    Atmosphere atmosphere_;
    std::optional<AtmosphereParameters> updated_parameters_{};

    GLBuffer buffer_;
    GLReloadableProgram program_flat_;
    GLReloadableProgram program_terrain_;
    GLReloadableComputeProgram program_shadow_;
    GLReloadableComputeProgram program_shadow_froxel_;
    GLSampler sampler_;

    using vec3 = glm::vec3;
    using vec2 = glm::vec2;
    struct TerrainUBO {
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
    } terrain_data_;
    GLBuffer terrain_ubo_;
};
