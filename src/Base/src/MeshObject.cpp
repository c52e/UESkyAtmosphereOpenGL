#include "MeshObject.h"

#include <imgui.h>
#include <glm/gtc/type_ptr.hpp>

#include "ShadowMap.h"
#include "GBuffer.h"

MeshObject::MeshObject(const char* name, const Mesh* mesh, const glm::mat4& model, bool cast_shadow)
    : name_(name)
    , mesh_(mesh)
    , model_(model)
    , cast_shadow_(cast_shadow) {}

void MeshObject::RenderToGBuffer(const glm::mat4& view_projection) const {
    if (visible_) {
        GBufferRenderer::Instance().Setup(model_, view_projection, material);
        mesh_->Draw();
    }
}

void MeshObject::RenderToShadowMap(const glm::mat4& light_view_projection) const {
    if (visible_ && cast_shadow_) {
        ShadowMapRenderer::Instance().Setup(model_, light_view_projection);
        mesh_->Draw();
    }
}

void MeshObject::DrawGui() {
    if (ImGui::TreeNode(name_.c_str())) {
        ImGui::Checkbox("Visible", &visible_);
        ImGui::SliderFloat("Metallic Factor", &material.metallic_factor, 0.0f, 1.0f);
        ImGui::SliderFloat("Roughness Factor", &material.roughness_factor, 0.0f, 1.0f);
        ImGui::SliderFloat("Ao Factor", &material.ao_factor, 0.0f, 1.0f);
        ImGui::ColorEdit3("Albedo Factor", glm::value_ptr(material.albedo_factor));
        ImGui::TreePop();
    }
}

Meshes::Meshes() {
    pyramid_ = std::make_unique<Mesh>(CreatePyramid());
    sphere_ = std::make_unique<Mesh>(CreateSphere());
    tyrannosaurus_rex_ = std::make_unique<Mesh>(ReadMeshFile("../data/models/TyrannosaurusRex.gltf"));
    wall_ = std::make_unique<Mesh>(ReadMeshFile("../data/models/wall.gltf"));
    fence_ = std::make_unique<Mesh>(ReadMeshFile("../data/models/fence.gltf"));
    cylinder_ = std::make_unique<Mesh>(ReadMeshFile("../data/models/cylinder.gltf"));
}
