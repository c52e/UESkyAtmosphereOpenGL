#pragma once

#include <vector>
#include <filesystem>

#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include "gl.hpp"

class GltfModel {
public:
	struct Material {
		glm::vec4 base_color_factor = glm::vec4(1.0f);
		glm::vec3 emissive_factor = glm::vec3(0.0f);
		float ao_factor = 1.0f;

		int base_color_index = -1;
		int normal_index = -1;
		int metallic_roughness_index = -1;
		int emissive_index = -1;
		int ao_index = -1;
	};

	struct Node {
		int mesh;
		std::vector<int> children;

		glm::vec3 translation{};
		glm::vec3 scale{ 1.0f };
		glm::quat rotation{};
		glm::mat4 matrix = glm::mat4(1.0f);

		glm::mat4 LocalMatrix() const {
			return glm::translate(glm::mat4(1.0f), translation) * glm::mat4(rotation) * glm::scale(glm::mat4(1.0f), scale) * matrix;
		}
	};

	struct Primitive {
		int mode;
		int material_index;
		int first_index;
		uint32_t count;
	};

	struct Mesh {
		std::vector<Primitive> primitives;
	};

	struct Vertex {
		glm::vec3 position;
		glm::vec3 normal;
		glm::vec2 uv;
		glm::vec3 tangent;
	};

	GltfModel(const std::filesystem::path& path);

	void UpdateMatrices(const glm::mat4& m);

	void RenderToGBuffer(const glm::mat4& view_projection) const;

	void RenderToShadowMap(const glm::mat4& light_view_projection) const;

private:
	GLuint GetTexture(int index) const {
		return index >= 0 ? images_[textures_[index]].id() : 0;
	}

	void UpdateMatrices(const glm::mat4& m, int i);

	std::vector<glm::mat4> matrices_;

	std::vector<Material> materials_;
	std::vector<GLTexture> images_;
	std::vector<int> textures_;
	std::vector<Node> nodes_;
	std::vector<Mesh> meshes_;
	std::vector<int> scene_nodes_;

	GLVertexArray vao_;
	GLBuffer vbo_;
	GLBuffer ebo_;
};
