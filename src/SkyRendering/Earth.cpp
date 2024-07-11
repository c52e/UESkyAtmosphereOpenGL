#include "Earth.h"

#include <array>

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Textures.h"
#include "Samplers.h"
#include "ScreenRectangle.h"
#include "Atmosphere.h"

using namespace glm;
#define saturate(x) clamp(x, 0.0f, 1.0f)
#include "../../shaders/SkyRendering/Terrain.glsl"

struct EarthBufferData {
    glm::mat4 view_projection;
    glm::mat4 inv_view_projection;
    glm::vec3 camera_position;
    float camera_earth_center_distance;
    glm::vec3 earth_center;
    float padding;
    glm::vec3 up_direction;
    float padding1;
    glm::vec3 light_dir;
    float padding2;
};

Earth::Earth() {
    buffer_.Create();
    glNamedBufferStorage(buffer_.id(), sizeof(EarthBufferData), NULL, GL_DYNAMIC_STORAGE_BIT);
    terrain_ubo_.Create();
    glNamedBufferStorage(terrain_ubo_.id(), sizeof(terrain_data_), NULL, GL_DYNAMIC_STORAGE_BIT);

    program_flat_ = []() {
        auto src = ReadWithPreprocessor("../shaders/SkyRendering/EarthRender.glsl");
        return GLProgram(kCommonVertexSrc, src.c_str());
    };
    program_terrain_ = []() {
        auto src = Replace(ReadWithPreprocessor("../shaders/SkyRendering/EarthRender.glsl"), "ENABLE_RAYMARCHED_TERRAIN", "1");
        return GLProgram(kCommonVertexSrc, src.c_str());
    };
    program_shadow_ = {
        "../shaders/SkyRendering/EarthRender.glsl",
        {{16, 8}, {8, 4}, {8, 8}, {16, 4}, {16, 16}, {32, 8}, {32, 16}},
        [](const std::string& src) { return Replace(src, "SHADOW_PROGRAM", "1"); }
    };
    program_shadow_froxel_ = {
        "../shaders/SkyRendering/EarthRender.glsl",
        {{16, 8}, {8, 4}, {8, 8}, {16, 4}, {16, 16}, {32, 8}, {32, 16}},
        [](const std::string& src) { return Replace(src, "SHADOW_FROXEL_PROGRAM", "1"); }
    };

    sampler_.Create();
    glSamplerParameteri(sampler_.id(), GL_TEXTURE_WRAP_S, GL_REPEAT);
    glSamplerParameteri(sampler_.id(), GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glSamplerParameteri(sampler_.id(), GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glSamplerParameteri(sampler_.id(), GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    float max_texture_max_anisotropy;
    glGetFloatv(GL_MAX_TEXTURE_MAX_ANISOTROPY, &max_texture_max_anisotropy);
    glSamplerParameterf(sampler_.id(), GL_TEXTURE_MAX_ANISOTROPY, max_texture_max_anisotropy);
}

void Earth::Update(const Camera& camera, glm::vec2 viewport, const glm::vec3& light_dir) {
    if (always_update_atmosphere
        || !updated_parameters_.has_value()
        || memcmp(&*updated_parameters_, &parameters, sizeof(parameters))) {
        atmosphere_.UpdateLuts(parameters);
        updated_parameters_.emplace(parameters);
    }

    EarthBufferData buffer;
    buffer.view_projection = camera.ViewProjection();
    buffer.inv_view_projection = glm::inverse(buffer.view_projection);
    buffer.camera_position = camera.position();
    buffer.earth_center = center();
    buffer.camera_earth_center_distance = glm::distance(buffer.camera_position, buffer.earth_center);
    buffer.up_direction = glm::normalize(buffer.camera_position - buffer.earth_center);
    buffer.light_dir = light_dir;

    glNamedBufferSubData(buffer_.id(), 0, sizeof(buffer), &buffer);

    auto pixel_tangent = camera.tangent_half_fovy() * 2.0f / static_cast<float>(viewport.y);
    auto base_lod_k = pixel_tangent / terrain.slope;
    terrain_data_.uAltitudeLodK = base_lod_k * glm::exp2(-(terrain.altitude_fineness * 20.0f - 10.0f));
    terrain_data_.uNormalLodK = base_lod_k * exp2(-(terrain.normal_fineness * 20.0f - 10.0f));
    terrain_data_.uMarchingThreshholdTangent = pixel_tangent * glm::exp2(-(terrain.marching_fineness * 20.0f - 10.0f));

    terrain_data_.uTerrainMaxAltitude = terrain.max_altitude;
    terrain_data_.uMaxVisibleDistance = terrain.max_visible_distance;
    terrain_data_.uNormalDeltaTangent = terrain.normal_delta_tangent;
    terrain_data_.uMarchingStepScale = 0.228f / terrain.slope;
    terrain_data_.uMaxSteps = terrain_data_.uMaxVisibleDistance / (0.5f * terrain_data_.uMarchingStepScale * terrain.max_altitude);
    terrain_data_.uTerrainFrequency = terrain.slope / terrain.max_altitude;
    terrain_data_.uDirtAlbedo = terrain.dirt_albedo;
    terrain_data_.uGrassAlbedo = terrain.grass_albedo;
    terrain_data_.uSnowAlbedo = terrain.snow_albedo;

    glNamedBufferSubData(terrain_ubo_.id(), 0, sizeof(terrain_data_), &terrain_data_);
}

glm::vec3 Earth::CheckCollision(const glm::vec3& target_pos, float zNear) {
    if (enable_raymarched_terrain) {
        auto miny = TerrainY(glm::vec2(target_pos.x, target_pos.z), terrain.physical_lod
            , terrain_data_.uTerrainFrequency, terrain_data_.uTerrainMaxAltitude
            , parameters.bottom_radius, center()) + zNear;
        if (target_pos.y < miny)
            return glm::vec3(target_pos.x, miny, target_pos.z);
    }
    else {
        auto radius = parameters.bottom_radius + zNear;
        if (glm::length(target_pos - center()) < radius)
            return glm::normalize(target_pos - center()) * (radius + 0.001f) + center();
    }
    return target_pos;
}

void Earth::RenderToGBuffer(GLuint depth_texture) {
    PERF_MARKER("Earth/Terrain Render");
    glBindBufferBase(GL_UNIFORM_BUFFER, 1, buffer_.id());
    glBindBufferBase(GL_UNIFORM_BUFFER, 2, terrain_ubo_.id());

    GLBindTextures({ depth_texture, Textures::Instance().earth_albedo() });
    GLBindSamplers({ 0u, sampler_.id() });

    glUseProgram((enable_raymarched_terrain ? program_terrain_ : program_flat_).id());
    glDepthFunc(GL_ALWAYS);
    ScreenRectangle::Instance().Draw();
    glDepthFunc(GL_LESS);
}

void Earth::RenderShadowBuffer(const GBuffer& gbuffer, const FroxelInfo& froxel_info) {
    if (!(enable_raymarched_terrain && enable_terrain_shadow))
        return;

    glBindBufferBase(GL_UNIFORM_BUFFER, 1, buffer_.id());
    glBindBufferBase(GL_UNIFORM_BUFFER, 2, terrain_ubo_.id());
    {
        PERF_MARKER("Terrain Shadow");
        GLBindTextures({ gbuffer.depth(), gbuffer.normal() });
        GLBindSamplers({ 0u, 0u });
        GLBindImageTextures({ gbuffer.pixel_visibility() });

        glUseProgram(program_shadow_.id());
        program_shadow_.Dispatch(gbuffer.dim());
    }
    {
        PERF_MARKER("Terrain Froxel Shadow");
        GLBindImageTextures({ froxel_info.froxel });
        glUseProgram(program_shadow_froxel_.id());
        glUniform1f(0, froxel_info.max_distance);
        program_shadow_froxel_.Dispatch(froxel_info.dim);
    }
    glMemoryBarrier(GL_TEXTURE_FETCH_BARRIER_BIT);
}

void Earth::DrawTerrainGui() {
    ImGui::Checkbox("Enable Terrain", &enable_raymarched_terrain);
    ImGui::SameLine();
    ImGui::Checkbox("Cast Shadow", &enable_terrain_shadow);

    ImGui::SliderFloat("Max Altitude", &terrain.max_altitude, 0.0f, 10.0f);
    ImGui::SliderFloat("Slope", &terrain.slope, 0.0f, 2.0f);
    ImGui::SliderFloat("Max Visible Distance", &terrain.max_visible_distance, 0.0f, 2000.0f);
    ImGui::SliderFloat("Normal Delta Tangent", &terrain.normal_delta_tangent, 0.0f, 0.005f, "%.5f"); 
    ImGui::SliderFloat("Altitude Fineness", &terrain.altitude_fineness, 0.0f, 1.0f);
    ImGui::SliderFloat("Normal Fineness", &terrain.normal_fineness, 0.0f, 1.0f);
    ImGui::SliderFloat("Marching Fineness", &terrain.marching_fineness, 0.0f, 1.0f);
    ImGui::SliderFloat("Physical Lod", &terrain.physical_lod, 0.0f, 12.0f);
    ImGui::ColorEdit3("Dirt Albedo", glm::value_ptr(terrain.dirt_albedo));
    ImGui::ColorEdit3("Grass Albedo", glm::value_ptr(terrain.grass_albedo));
    ImGui::ColorEdit3("Snow Albedo", glm::value_ptr(terrain.snow_albedo));
}

glm::mat4 Earth::moon_model() const {
    auto model = glm::identity<glm::mat4>();
    model = glm::translate(model, center());
    model = glm::rotate(model, glm::radians(moon_status.direction_phi), glm::vec3(0, 1, 0));
    model = glm::rotate(model, glm::radians(moon_status.direction_theta), glm::vec3(0, 0, -1));
    model = glm::translate(model,
        glm::vec3(0, moon_status.distance + moon_status.radius + parameters.bottom_radius, 0));
    model = glm::scale(model, glm::vec3(moon_status.radius));
    model = glm::rotate(model, glm::radians(90.f), glm::vec3(0, 0, 1));
    return model;
}
