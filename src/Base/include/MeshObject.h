#pragma once

#include <string>
#include <memory>

#include <glad/glad.h>
#include <glm/glm.hpp>

#include "Mesh.h"
#include "Singleton.h"

struct Material {
    glm::vec3 albedo_factor{ 1, 1, 1};
    GLuint albedo_texture = 0;
    float metallic_factor = 0.0f;
    float roughness_factor = 1.0f;
    float ao_factor = 1.0f;
    GLuint normal_texture = 0;
    GLuint metallic_roughness_texture = 0;
    GLuint ao_texture = 0;
};

class MeshObject {
public:
    MeshObject(const char* name, const Mesh* mesh, const glm::mat4& model, bool cast_shadow);

    void RenderToGBuffer(const glm::mat4& view_projection) const;

    void RenderToShadowMap(const glm::mat4& light_view_projection) const;

    const glm::mat4& model() const { return model_; }
    void set_model(const glm::mat4& model) { model_ = model; }

    void DrawGui();

    bool visible_ = true;
    Material material;

protected:
    std::string name_;
    const Mesh* mesh_;
    glm::mat4 model_;
    bool cast_shadow_;
};

class Meshes :public Singleton<Meshes> {
public:
	friend Singleton<Meshes>;

    Mesh* pyramid() const {
		return pyramid_.get();
	}

    Mesh* sphere() const {
        return sphere_.get();
    }

    Mesh* tyrannosaurus_rex() const {
        return tyrannosaurus_rex_.get();
    }

    Mesh* wall() const {
        return wall_.get();
    }

    Mesh* fence() const {
        return fence_.get();
    }

    Mesh* cylinder() const {
        return cylinder_.get();
    }

private:
    Meshes();

    std::unique_ptr<Mesh> pyramid_;
    std::unique_ptr<Mesh> sphere_;
    std::unique_ptr<Mesh> tyrannosaurus_rex_;
    std::unique_ptr<Mesh> wall_;
    std::unique_ptr<Mesh> fence_;
    std::unique_ptr<Mesh> cylinder_;
};
