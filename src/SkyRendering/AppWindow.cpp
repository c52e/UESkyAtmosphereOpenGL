#include "AppWindow.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>

#include <magic_enum.hpp>
#include <imgui.h>
#include <rapidjson/error/en.h>

#include "Textures.h"
#include "Samplers.h"
#include "ImGuiExt.h"
#include "ScreenRectangle.h"
#include "Profiler.h"

AppWindow::AppWindow(const char* config_path, int width, int height)
    : GLWindow((std::string("SkyRendering (") + config_path + ")").c_str(), width, height, false) {

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    auto aspect = static_cast<float>(width) / height;
    camera_ = Camera(glm::vec3(0.0f, 1.0f, -1.0f), 180.0f, 0.0f, aspect);
    camera_.zNear = 3e-1f;
    HandleReshapeEvent(width, height);
    shadow_map_ = std::make_unique<ShadowMap>(2048, 2048);

    Init(config_path);
}

void AppWindow::Init(const char* config_path) {
    glClipControl(GL_LOWER_LEFT, GL_ZERO_TO_ONE);

    std::ifstream fin(config_path);
    if (fin) {
        using namespace rapidjson;
        auto str = std::string(std::istreambuf_iterator<char>{fin}, {});
        Document d;
        if (d.Parse<kParseCommentsFlag | kParseTrailingCommasFlag>(str.c_str()).HasParseError()) {
            std::ostringstream str;
            str << "Failed to parse \"" << config_path << "\" (" << "offset " << d.GetErrorOffset() << "): " << GetParseError_En(d.GetParseError());
            throw std::runtime_error(str.str());
        }
        Deserialize(*this, d);
        LOG_INFO("\"{}\" loaded", config_path);
    }
    else {
        SaveConfig(config_path);
        LOG_WARN("\"{}\" not found. Default config created", config_path);
    }
    strcpy_s(config_path_, config_path);

    SetFullScreen(full_screen_);
    Samplers::SetAnisotropyEnable(anisotropy_enable_);
    glfwSwapInterval(vsync_enable_ ? 1 : 0);
    atmosphere_renderer_ = std::make_unique<AtmosphereRenderer>(atmosphere_render_init_parameters_);

    mesh_objects_.clear();
    {
        auto object = std::make_unique<MeshObject>("Moon", Meshes::Instance().sphere(), earth_.moon_model(), false);
        object->material.albedo_texture = Textures::Instance().moon_albedo();
        object->material.normal_texture = Textures::Instance().moon_normal();
        moon_ = object.get();
        mesh_objects_.push_back(std::move(object));
    }
    {
        auto model = glm::identity<glm::mat4>();
        model = glm::translate(model, glm::vec3(0, 1, 1));
        auto object = std::make_unique<MeshObject>("Tiny Earth", Meshes::Instance().sphere(), model, true);
        object->material.albedo_texture = Textures::Instance().earth_albedo();
        mesh_objects_.push_back(std::move(object));
    }
    {
        auto model = glm::identity<glm::mat4>();
        model = glm::translate(model, glm::vec3(0, 1, 3));
        model = glm::scale(model, glm::vec3(0.4f));
        auto object = std::make_unique<MeshObject>("Tiny Moon", Meshes::Instance().sphere(), model, true);
        object->material.albedo_texture = Textures::Instance().moon_albedo();
        object->material.normal_texture = Textures::Instance().moon_normal();
        mesh_objects_.push_back(std::move(object));
    }
    {
        auto model = glm::identity<glm::mat4>();
        model = glm::translate(model, glm::vec3(0, 0, -1));
        model = glm::rotate(model, glm::pi<float>() / 2, glm::vec3(0, -1, 0));
        model = glm::rotate(model, glm::pi<float>() / 2, glm::vec3(-1, 0, 0));
        model = glm::scale(model, glm::vec3(0.3e-3, 0.3e-3, 0.3e-3));
        auto object = std::make_unique<MeshObject>("Tyrannosaurus", Meshes::Instance().tyrannosaurus_rex(), model, true);
        object->material.albedo_factor = { 0.42, 0.42, 0.42 };
        mesh_objects_.push_back(std::move(object));
    }
    {
        auto model = glm::identity<glm::mat4>();
        model = glm::translate(model, glm::vec3(-3, 0, -1.5));
        model = glm::rotate(model, glm::pi<float>() / 2, glm::vec3(0, -1, 0));
        model = glm::scale(model, glm::vec3(1e-1, 1e-1, 1e-1));
        auto object = std::make_unique<MeshObject>("Wall", Meshes::Instance().wall(), model, true);
        object->material.albedo_factor = { 0.42, 0.42, 0.42 };
        mesh_objects_.push_back(std::move(object));
    }
    {
        auto model = glm::identity<glm::mat4>();
        model = glm::translate(model, glm::vec3(-2, -0.1, -2));
        model = glm::rotate(model, glm::pi<float>() / 2, glm::vec3(-1, 0, 0));
        model = glm::scale(model, glm::vec3(1e-2, 1e-2, 1e-2));
        auto object = std::make_unique<MeshObject>("Fence", Meshes::Instance().fence(), model, true);
        object->material.albedo_factor = { 0.42, 0.42, 0.42 };
        mesh_objects_.push_back(std::move(object));
    }
    {
        auto model = glm::identity<glm::mat4>();
        model = glm::translate(model, glm::vec3(-2.2, 0.8, -2.0));
        model = glm::rotate(model, glm::pi<float>() / 2, glm::vec3(0, 1, 0));
        model = glm::rotate(model, glm::pi<float>() / 2, glm::vec3(0, 0, 1));
        model = glm::scale(model, glm::vec3(1e-2, 1e-2, 1e-2));
        auto object = std::make_unique<MeshObject>("Fence2", Meshes::Instance().fence(), model, true);
        object->material.albedo_factor = { 0.42, 0.42, 0.42 };
        mesh_objects_.push_back(std::move(object));
    }
    {
        auto model = glm::identity<glm::mat4>();
        model = glm::translate(model, glm::vec3(2, 0, 0));
        model = glm::rotate(model, glm::pi<float>() / 2, glm::vec3(-1, 0, 0));
        model = glm::scale(model, glm::vec3(1e-1, 1e-1, 1e1));
        auto object = std::make_unique<MeshObject>("Cylinder", Meshes::Instance().cylinder(), model, true);
        object->material.albedo_factor = { 0.42, 0.42, 0.42 };
        mesh_objects_.push_back(std::move(object));
    }
}

void AppWindow::HandleDropEvent(int count, const char** paths) {
    gltf_models_.emplace_back(std::make_unique<GltfModel>(paths[0]));
}

bool AppWindow::SaveConfig(const char* path) {
    rapidjson::StringBuffer sb;
    rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(sb);
    Serialize(*this, writer);
    std::ofstream fout(path);
    if (!fout)
        return false;
    fout << sb.GetString();
    return true;
}

void AppWindow::HandleDisplayEvent() {
    auto [width, height] = GetWindowSize();
    earth_.Update(camera_, { width, height }, atmosphere_renderer_->sun_direction());
    camera_.set_position(earth_.CheckCollision(camera_.position(), camera_.zNear));
    volumetric_cloud_.Update(camera_, earth_, *atmosphere_renderer_);
    moon_->set_model(earth_.moon_model());


    auto model = glm::identity<glm::mat4>();
    model = glm::translate(model, glm::vec3(0, 1, -4));
    for (auto& gltf: gltf_models_)
        gltf->UpdateMatrices(model);
    
    Render();
    ProcessInput();
}

void AppWindow::Render() {
    // https://www.reddit.com/r/programming/comments/4tlr2u/john_carmack_on_inlined_code/
    // http://number-none.com/blow/john_carmack_on_inlined_code.html
    PERF_MARKER("Render")

    glEnable(GL_FRAMEBUFFER_SRGB);
    constexpr float kShadowRegionHalfWidth = 4.0f;
    constexpr float kShadowRegionNear = 0.0f;
    constexpr float kShadowRegionFar = 5e1f;

    auto light_view_projection_matrix = ComputeLightMatrix(5.f,
        atmosphere_render_parameters_.sun_direction_theta,
        atmosphere_render_parameters_.sun_direction_phi,
        -kShadowRegionHalfWidth, kShadowRegionHalfWidth,
        -kShadowRegionHalfWidth, kShadowRegionHalfWidth,
        kShadowRegionNear, kShadowRegionFar);

    {
        PERF_MARKER("RenderShadowMap");
        glEnable(GL_DEPTH_TEST);

        glCullFace(GL_FRONT);
        glEnable(GL_CULL_FACE);
        shadow_map_->ClearBindViewport();
        for (const auto& mesh_object : mesh_objects_)
            mesh_object->RenderToShadowMap(light_view_projection_matrix);
        for (auto& gltf : gltf_models_)
            gltf->RenderToShadowMap(light_view_projection_matrix);
        glDisable(GL_DEPTH_TEST);
    }

    auto [width, height] = GetWindowSize();
    glViewport(0, 0, width, height);
    {
        PERF_MARKER("RenderGBuffer");
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_GREATER);

        glBindFramebuffer(GL_FRAMEBUFFER, gbuffer_->id());
        Clear(*gbuffer_);
        glCullFace(GL_BACK);
        for (const auto& mesh_object : mesh_objects_)
            mesh_object->RenderToGBuffer(camera_.ViewProjection());
        for (auto& gltf : gltf_models_)
            gltf->RenderToGBuffer(camera_.ViewProjection());
        earth_.RenderToGBuffer(gbuffer_->depth());
        glDepthFunc(GL_LESS);
        glDisable(GL_DEPTH_TEST);
    }
    {
        volumetric_cloud_.RenderShadow(*gbuffer_);
        auto froxel = volumetric_cloud_.GetShadowFroxel();
        earth_.RenderShadowBuffer(*gbuffer_, { froxel.shadow_froxel, froxel.dim, froxel.max_distance});
    }
    hdrbuffer_->BindHdrFramebuffer();
    {
        PERF_MARKER("RenderViewport");
        auto sun_angular_radius = glm::radians(earth_.parameters.sun_angular_radius);
        atmosphere_render_parameters_.pcss_size_k = sun_angular_radius * (kShadowRegionFar - kShadowRegionNear) / kShadowRegionHalfWidth;
        atmosphere_render_parameters_.blocker_kernel_size_k = 2.0f * atmosphere_render_parameters_.pcss_size_k;
        atmosphere_render_parameters_.view_projection = camera_.ViewProjection();
        atmosphere_render_parameters_.camera_position = camera_.position();
        atmosphere_render_parameters_.light_view_projection = light_view_projection_matrix;
        atmosphere_render_parameters_.shadow_map_texture = shadow_map_->depth_texture();
        atmosphere_render_parameters_.gbuffer = gbuffer_.get();
        atmosphere_render_parameters_.hdrbuffer = hdrbuffer_.get();

        atmosphere_renderer_->Render(earth_, volumetric_cloud_, atmosphere_render_parameters_);
    }
    {
        volumetric_cloud_.Render(hdrbuffer_->hdr_texture(), gbuffer_->depth());
    }

    hdrbuffer_->DoPostProcessAndBindSdrFramebuffer(post_process_parameters_);
    smaa_->DoSMAA(hdrbuffer_->sdr_texture());
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    TextureVisualizer::Instance().VisualizeTexture(smaa_->output_tex());
    glDisable(GL_FRAMEBUFFER_SRGB);
}

void AppWindow::ProcessInput() {
    if (ImGui::GetIO().WantCaptureKeyboard || ImGui::GetIO().WantCaptureMouse)
        return;

    static float previous_frame_time_ = 0.0f;
    float frame_time = static_cast<float>(glfwGetTime());
    float dt = frame_time - previous_frame_time_;
    previous_frame_time_ = frame_time;

    float dForward = 0.0f;
    float dRight = 0.0f;
    float dUp = 0.0f;
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        dForward += camera_speed_ * dt;
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        dForward -= camera_speed_ * dt;
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        dRight += camera_speed_ * dt;
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        dRight -= camera_speed_ * dt;
    if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS)
        dUp += camera_speed_ * dt;
    if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
        dUp -= camera_speed_ * dt;
    auto move_vector = camera_.front() * dForward + camera_.right() * dRight + kWorldUp * dUp;
    auto new_position = camera_.position() + move_vector;
    camera_.set_position(new_position);

    float rotate_speed = 90.0f;
    float dPitch = 0.0f;
    float dYaw = 0.0f;
    if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS)
        dPitch += rotate_speed * dt;
    if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS)
        dPitch -= rotate_speed * dt;
    if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS)
        dYaw += rotate_speed * dt;
    if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS)
        dYaw -= rotate_speed * dt;

    static bool previoud_left_button_down = false;
    static double previous_mouse_x_ = 0;
    static double previous_mouse_y_ = 0;
    bool left_button_down = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS;
    if (!previoud_left_button_down && left_button_down) {
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
        previous_mouse_y_ = mouse_y_;
        previous_mouse_x_ = mouse_x_;
    }
    else if (previoud_left_button_down && !left_button_down) {
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    }
    previoud_left_button_down = left_button_down;

    if (left_button_down) {
        float mouse_rotate_speed = 10.f;
        dPitch += mouse_rotate_speed * static_cast<float>(previous_mouse_y_ - mouse_y_) * dt;
        dYaw += mouse_rotate_speed * static_cast<float>(mouse_x_ - previous_mouse_x_) * dt;
        previous_mouse_y_ = mouse_y_;
        previous_mouse_x_ = mouse_x_;
    }

    camera_.Rotate(dPitch, dYaw);
}

static void SliderFloat(const char* name, float* ptr, float min_v, float max_v) {
    ImGui::SliderFloat(name, ptr, min_v, max_v);
}

static void SliderFloatLogarithmic(const char* name, float* ptr, float min_v, float max_v, const char* format = "%.3f") {
    ImGui::SliderFloat(name, ptr, min_v, max_v, format, ImGuiSliderFlags_Logarithmic);
}

static void ColorEdit(const char* name, glm::vec3& color) {
    ImGui::ColorEdit3(name, glm::value_ptr(color));
}

void AppWindow::HandleDrawGuiEvent() {
    if (!draw_gui_enable_)
        return;
    ImGui::Begin("Main");
    ImGui::Text("%.3f ms (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);

    if (ImGui::Button("Reload Shader")) {
        try {
            GLReloadableProgram::ReloadAll();
        }
        catch (std::exception& e) {
            LOG_ERROR("{}", e.what());
        }
    }
    ImGui::SameLine();
    {
        static std::vector<std::string> config_paths;
        static std::vector<const char*> config_paths_cstr;
        static int index = 0;
        if (ImGui::Button("Load State..")) {
            ImGui::OpenPopup("Load State");
            namespace fs = std::filesystem;
            config_paths.clear();
            config_paths_cstr.clear();
            for (auto const& entry : fs::recursive_directory_iterator(fs::current_path())) {
                if (fs::is_regular_file(entry) && entry.path().extension() == ".json") {
                    config_paths.emplace_back(fs::relative(entry.path(), fs::current_path()).string());
                    if (config_paths.back() == config_path_)
                        index = static_cast<int>(config_paths.size()) - 1;
                }
            }
            for (const auto& path : config_paths)
                config_paths_cstr.emplace_back(path.c_str());
        }
        if (ImGui::BeginPopupModal("Load State", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) {
            if (!config_paths_cstr.empty()) {
                if (index >= static_cast<int>(config_paths_cstr.size()))
                    index = static_cast<int>(config_paths_cstr.size()) - 1;
                ImGui::ListBox("Configs", &index, config_paths_cstr.data(), static_cast<int>(config_paths_cstr.size()));
                if (ImGui::Button("OK", ImVec2(120, 0))) {
                    try {
                        Init(config_paths_cstr[index]);
                    }
                    catch (std::exception& e) {
                        Error(std::string("Load state failed: ") + e.what());
                    }
                    ImGui::CloseCurrentPopup();
                }
            }
            ImGui::SameLine();
            if (ImGui::Button("Cancel", ImVec2(120, 0))) {
                ImGui::CloseCurrentPopup();
            }
            ImGui::EndPopup();
        }
    }
    
    ImGui::SameLine();
    if (ImGui::Button("Save State..")) {
        ImGui::OpenPopup("Save State");
    }
    if (ImGui::BeginPopupModal("Save State", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) {
        ImGui::InputText("path", config_path_, std::size(config_path_));
        if (ImGui::Button("OK", ImVec2(120, 0))) {
            if (!SaveConfig(config_path_))
                Error("Save failed");
            ImGui::CloseCurrentPopup();
        }
        ImGui::SameLine();
        if (ImGui::Button("Cancel", ImVec2(120, 0))) {
            ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
    }
    
    ImGui::SameLine();
    static ImageFormat format = ImageFormat::PNG;
    if (ImGui::Button("ScreenShot") || need_screenshot_) {
        time_t rawtime{};
        std::time(&rawtime);
        std::tm timeinfo{};
        localtime_s(&timeinfo, &rawtime);
        char buffer[84];
        std::strftime(buffer, std::size(buffer), "%Y-%m-%d %H-%M-%S", &timeinfo);

        namespace fs = std::filesystem;
        std::string folder = "screenshots/";
        if (!fs::exists(folder))
            fs::create_directory(folder);
        ScreenShot((folder + buffer).c_str(), format);
        need_screenshot_ = false;
    }
    ImGui::SameLine();
    ImGui::PushItemWidth(60);
    ImGui::EnumSelect("Type", &format);
    ImGui::PopItemWidth();

    ImGui::SameLine();
    static int value_index = 11;
    static int interval = 50;
    if (ImGui::Button("ScreenShotRayleighMie")) {
        value_index = 0;
    }
    if (value_index < 11) {
        auto r_index = value_index / 3;
        auto m_index = value_index % 3;
        auto o_index = 1;
        if (value_index >= 9) {
            r_index = 1;
            m_index = 1;
            o_index = value_index == 9 ? 0 : 2;
        }
        earth_.parameters.rayleigh_scattering_scale = std::array{ 0.0331f * 0.33f, 0.0331f, 0.0331f * 3.0f }[r_index];
        earth_.parameters.mie_scattering_scale = std::array{ 0.003996f * 0.1f, 0.003996f, 0.003996f * 10.0f }[m_index];
        earth_.parameters.ozone_absorption_scale = std::array{ 0.0f, 0.001881f,  0.001881f * 5.0f}[o_index];
        if (--interval == 0) {
            interval = 50;
            char buffer[256];
            if (value_index < 9) {
                snprintf(buffer, std::size(buffer), "screenshots/R%dM%d", r_index, m_index);
            }
            else {
                snprintf(buffer, std::size(buffer), "screenshots/O%d", o_index);
            }
            ScreenShot(buffer, ImageFormat::JPG);
            value_index++;
        }
    }

    if (ImGui::Button("Hide GUI"))
        draw_gui_enable_ = false;
    ImGui::SameLine();
    ImGui::Text("(Press G to reshow)");
    ImGui::SameLine();
    ImGui::Checkbox("Show Debug Textures", &draw_debug_textures_enable_);
    ImGui::SameLine();
    ImGui::Checkbox("Show Help", &draw_help_enable_);

    bool previous_full_screen = full_screen_;
    bool previous_anisotropy_enable = anisotropy_enable_;
    bool previous_vsync_enable = vsync_enable_;
    auto previous_atmosphere_render_init_parameters = atmosphere_render_init_parameters_;
    auto previous_smaa_option = smaa_option_;
    auto previous_volumetric_cloud_full_resolution = volumetric_cloud_full_resolution_;

    ImGui::Checkbox("Full Screen", &full_screen_);
    ImGui::SameLine();
    ImGui::Checkbox("VSync", &vsync_enable_);

    if (ImGui::TreeNode("Render")) {
        ImGui::Checkbox("PCSS", &atmosphere_render_init_parameters_.pcss_enable);
        ImGui::SameLine();
        ImGui::Checkbox("Volumetric Light", &atmosphere_render_init_parameters_.volumetric_light_enable);
        ImGui::SameLine();
        ImGui::Checkbox("Moon Shadow", &atmosphere_render_init_parameters_.moon_shadow_enable);
        ImGui::SameLine();
        ImGui::Checkbox("Anisotropy Filtering", &anisotropy_enable_);
        ImGui::EnumSelect("SMAA", &smaa_option_);
        ImGui::Separator();

        SliderFloat("Transmittance Steps", &earth_.parameters.transmittance_steps, 0, 100.0);
        SliderFloat("Multiscattering Steps", &earth_.parameters.multiscattering_steps, 0, 100.0);
        SliderFloat("Multiscattering Mask", &earth_.parameters.multiscattering_mask, 0, 1.0);
        ImGui::Separator();

        ImGui::Checkbox("Use Sky View LUT", &atmosphere_render_init_parameters_.use_sky_view_lut);
        ImGui::Checkbox("Sky View LUT Dither Sample Point Enable", &atmosphere_render_init_parameters_.sky_view_lut_dither_sample_point_enable);
        SliderFloat("Sky View LUT Steps", &atmosphere_render_parameters_.sky_view_lut_steps, 0, 100.0);
        ImGui::Separator();

        ImGui::Checkbox("Use Aerial Perspective LUT", &atmosphere_render_init_parameters_.use_aerial_perspective_lut);
        ImGui::Checkbox("Aerial Perspective LUT Dither Sample Point Enable", &atmosphere_render_init_parameters_.aerial_perspective_lut_dither_sample_point_enable);
        SliderFloat("Aerial Perspective LUT Steps", &atmosphere_render_parameters_.aerial_perspective_lut_steps, 0, 100.0);
        SliderFloatLogarithmic("Aerial Perspective LUT Max Distance", &atmosphere_render_parameters_.aerial_perspective_lut_max_distance, 0, 7000.0, "%.0f");
        ImGui::SliderInt("Aerial Perspective LUT Depth", &atmosphere_render_init_parameters_.aerial_perspective_lut_depth, 1, 512);
        ImGui::Separator();

        ImGui::Checkbox("Raymarching Dither Sample Point Enable", &atmosphere_render_init_parameters_.raymarching_dither_sample_point_enable);
        SliderFloat("Raymarching Steps", &atmosphere_render_parameters_.raymarching_steps, 0, 100.0);
        ImGui::Checkbox("Use Scattering LUT", &atmosphere_render_init_parameters_.use_scattering_lut);
        ImGui::SliderInt("Scattering Order", &earth_.parameters.num_scattering_orders, 0, 10);
        ImGui::TreePop();
    }

    if (ImGui::TreeNode("Camera")) {
        SliderFloatLogarithmic("zNear", &camera_.zNear, 1e-3f, 1e1f, "%.3f");
        SliderFloatLogarithmic("fovy", &camera_.fovy, 0.0f, 180.0, "%.1f");
        SliderFloatLogarithmic("Move Speed", &camera_speed_, 0.f, 1e5f, "%.1f");
        ImGui::TreePop();
    }

    if (ImGui::TreeNode("Sun & Moon & Star")) {
        ColorEdit("Solar Illuminance", earth_.parameters.solar_illuminance);
        SliderFloat("Sun Angular Radius", &earth_.parameters.sun_angular_radius, 0, 30);
        SliderFloat("Sun Direction Theta", &atmosphere_render_parameters_.sun_direction_theta, 0, 180);
        SliderFloat("Sun Direction Phi", &atmosphere_render_parameters_.sun_direction_phi, 0, 360);
        ImGui::Separator();

        SliderFloat("Moon Direction Theta", &earth_.moon_status.direction_theta, 0, 180);
        SliderFloat("Moon Direction Phi", &earth_.moon_status.direction_phi, 0, 360);
        SliderFloatLogarithmic("Moon Distance From Earth", &earth_.moon_status.distance, 0, 1e6f, "%.0f");
        SliderFloatLogarithmic("Moon Radius", &earth_.moon_status.radius, 0, 1e4f, "%.0f");
        ImGui::Separator();

        SliderFloatLogarithmic("Star Luminance Scale", &atmosphere_render_parameters_.star_luminance_scale, 0, 1.0);
        ImGui::TreePop();
    }

    if (ImGui::TreeNode("Earth")) {
        ImGui::Checkbox("Always Update Atmosphere", &earth_.always_update_atmosphere);
        ColorEdit("Ground Albedo", earth_.parameters.ground_albedo);
        SliderFloatLogarithmic("Bottom Radius", &earth_.parameters.bottom_radius, 0, 6360.0, "%.1f");
        SliderFloatLogarithmic("Thickness", &earth_.parameters.thickness, 0, 1000.0, "%.1f");
        ImGui::Separator();

        SliderFloatLogarithmic("Rayleigh Exponential Distribution", &earth_.parameters.rayleigh_exponential_distribution, 0, 1000.0, "%.1f");
        SliderFloatLogarithmic("Rayleigh Scattering Scale", &earth_.parameters.rayleigh_scattering_scale, 0, 1.0f, "%.4f");
        ColorEdit("Rayleigh Scattering", earth_.parameters.rayleigh_scattering);
        ImGui::Separator();

        SliderFloatLogarithmic("Mie Exponential Distribution", &earth_.parameters.mie_exponential_distribution, 0, 1000.0, "%.1f");
        SliderFloat("Mie Phase G", &earth_.parameters.mie_phase_g, -1.0, 1.0);
        SliderFloatLogarithmic("Mie Scattering Scale", &earth_.parameters.mie_scattering_scale, 0, 1.0f, "%.4f");
        ColorEdit("Mie Scattering", earth_.parameters.mie_scattering);
        SliderFloatLogarithmic("Mie Absorption Scale", &earth_.parameters.mie_absorption_scale, 0, 1.0f, "%.4f");
        ColorEdit("Mie Absorption", earth_.parameters.mie_absorption);
        ImGui::Separator();
        
        SliderFloatLogarithmic("Ozone Center Altitude", &earth_.parameters.ozone_center_altitude, 0, 1000.0, "%.1f");
        SliderFloatLogarithmic("Ozone Width", &earth_.parameters.ozone_width, 0, 1000.0, "%.1f");
        SliderFloatLogarithmic("Ozone Absorption Scale", &earth_.parameters.ozone_absorption_scale, 0, 1.0f, "%.4f");
        ColorEdit("Ozone Absorption", earth_.parameters.ozone_absorption);
        ImGui::Separator();

        ImGui::TreePop();
    }
    if (ImGui::TreeNode("Terrain")) {
        earth_.DrawTerrainGui();
        ImGui::TreePop();
    }

    if (ImGui::TreeNode("Post Process")) {
        ImGui::EnumSelect("Tone Mapping", &post_process_parameters_.tone_mapping);
        ImGui::Checkbox("Dither Color", &post_process_parameters_.dither_color_enable);
        SliderFloat("Bloom Min Luminance", &post_process_parameters_.bloom_min_luminance, 0, 2.0f);
        SliderFloat("Bloom Max Clamp Luminance", &post_process_parameters_.bloom_max_delta_luminance, 0, 50.0f);
        SliderFloat("Bloom Filter Width", &post_process_parameters_.bloom_filter_width, 0, 0.1f);
        SliderFloat("Bloom Intensity", &post_process_parameters_.bloom_intensity, 0, 1.0f);
        SliderFloat("Exposure", &post_process_parameters_.exposure, 0, 100.0f);
        ImGui::TreePop();
    }

    if (ImGui::TreeNode("VolumetricCloud")) {
        ImGui::Checkbox("Full Resolution", &volumetric_cloud_full_resolution_);
        volumetric_cloud_.DrawGUI();
        ImGui::TreePop();
    }

    ImGui::End();

    if (draw_debug_textures_enable_) {
        auto show_image = [](GLuint tex, const char* name) {
            if (ImGui::TreeNode(name)) {
                ImGui::Image((void*)(intptr_t)tex, ImVec2(512, 512), ImVec2(0, 1), ImVec2(1, 0));
                ImGui::TreePop();
            }
        };
        ImGui::Begin("Textures");
        show_image(earth_.atmosphere().transmittance_texture(), "Transmittance LUT");
        show_image(earth_.atmosphere().multiscattering_texture(), "Multiscattering LUT");
        show_image(gbuffer_->albedo(), "Albedo");
        show_image(gbuffer_->normal(), "Normal");
        show_image(gbuffer_->orm(), "ORM");
        show_image(gbuffer_->depth(), "Depth");
        show_image(shadow_map_->depth_texture(), "Shadow Map Depth");
        show_image(Textures::Instance().env_brdf_lut(), "Env BRDF LUT");
        show_image(gbuffer_->pixel_visibility(), "Pixel Visibility");
        ImGui::End();
    }

    ImGui::Begin("Models");
    for (const auto& object : mesh_objects_) {
        object->DrawGui();
    }
    ImGui::End();

    ImGui::Begin("Compute Program");
    GLReloadableComputeProgram::DrawGUIAll();
    ImGui::End();

    if (draw_help_enable_) {
        ImGui::Begin("Help");
        ImGui::BulletText("The main window may be hidden by others. You can drag them away to show the main window and change parameters.");
        ImGui::BulletText("WASDQE to move camera (You can change speed in main window).");
        ImGui::BulletText("Hold right mouse button to rotate camera.");
        ImGui::End();
    }

    ImGui::Begin("Profile");
    Profiler::Instance().DrawImGui();
    ImGui::End();

    if (full_screen_ != previous_full_screen)
        SetFullScreen(full_screen_);

    if (vsync_enable_ != previous_vsync_enable)
        glfwSwapInterval(vsync_enable_ ? 1 : 0);

    if (anisotropy_enable_ != previous_anisotropy_enable)
        Samplers::SetAnisotropyEnable(anisotropy_enable_);

    if (atmosphere_render_init_parameters_ != previous_atmosphere_render_init_parameters)
        atmosphere_renderer_ = std::make_unique<AtmosphereRenderer>(atmosphere_render_init_parameters_);

    auto [width, height] = GetWindowSize();
    if (smaa_option_ != previous_smaa_option)
        smaa_ = std::make_unique<SMAA>(width, height, smaa_option_);

    if (volumetric_cloud_full_resolution_ != previous_volumetric_cloud_full_resolution)
        volumetric_cloud_.SetViewport(width, height, volumetric_cloud_full_resolution_);
}

void AppWindow::HandleReshapeEvent(int viewport_width, int viewport_height) {
    if (viewport_width > 0 && viewport_height > 0) {
        volumetric_cloud_.SetViewport(viewport_width, viewport_height, volumetric_cloud_full_resolution_);
        gbuffer_ = std::make_unique<GBuffer>(viewport_width, viewport_height);
        hdrbuffer_ = std::make_unique<HDRBuffer>(viewport_width, viewport_height);
        camera_.set_aspect(static_cast<float>(viewport_width) / viewport_height);
        smaa_ = std::make_unique<SMAA>(viewport_width, viewport_height, smaa_option_);
    }
}

void AppWindow::HandleKeyboardEvent(int key) {
    switch (key) {
    case GLFW_KEY_ESCAPE:
        Close();
        break;
    case GLFW_KEY_G:
        draw_gui_enable_ = !draw_gui_enable_;
        break;
    case GLFW_KEY_SPACE:
        need_screenshot_ = true;
        break;
    }
}

void AppWindow::HandleMouseEvent(double mouse_x, double mouse_y) {
    mouse_x_ = mouse_x;
    mouse_y_ = mouse_y;
}
