#include "GltfModel.h"

#include <unordered_set>

#include <glm/gtc/type_ptr.hpp>
#include <tiny_gltf.h>

#include "log.h"
#include "ImageLoader.h"
#include "GBuffer.h"
#include "ShadowMap.h"

namespace fs = std::filesystem;

using Node = GltfModel::Node;
using Primitive = GltfModel::Primitive;
using Vertex = GltfModel::Vertex;

std::vector<GltfModel::Material> LoadMaterials(const tinygltf::Model& input) {
    std::vector<GltfModel::Material> materials;
    materials.reserve(input.materials.size());
    for (const auto& material : input.materials) {
        GltfModel::Material m;
        if (auto itr = material.values.find("baseColorFactor"); itr != material.values.end()) {
            m.base_color_factor = glm::make_vec4(itr->second.ColorFactor().data());
        }
        if (auto itr = material.values.find("baseColorTexture"); itr != material.values.end()) {
            m.base_color_index = itr->second.TextureIndex();
        }

        m.normal_index = material.normalTexture.index;
        m.metallic_roughness_index = material.pbrMetallicRoughness.metallicRoughnessTexture.index;
        m.emissive_index = material.emissiveTexture.index;
        for (int j = 0; j < 3; ++j)
            m.emissive_factor[j] = material.emissiveFactor[j];
        m.ao_index = material.occlusionTexture.index;
        m.ao_factor = material.occlusionTexture.strength;

        materials.push_back(m);
    }

    return materials;
}

std::vector<GLTexture> LoadImages(const tinygltf::Model& input, const std::unordered_set<int>& srgb_indices) {
    std::vector<GLTexture> images;
    images.reserve(input.images.size());
    for (const auto& image : input.images) {
        ASSERT(image.bits == 8);

        auto width = image.width;
        auto height = image.height;
        auto component = image.component;
        auto pbuffer = &image.image;
        std::vector<uint8_t> rgba_buffer;
        if (component == 3) {
            rgba_buffer.resize(width * height * 4);
            auto src = image.image.data();
            auto dst = rgba_buffer.data();
            for (size_t i = 0; i < width * height; ++i) {
                *dst++ = *src++;
                *dst++ = *src++;
                *dst++ = *src++;
                *dst++ = 255;
            }
            component = 4;
            pbuffer = &rgba_buffer;
        }
        ASSERT(component == 4);

        GLTexture texture;
        texture.Create(GL_TEXTURE_2D);
        auto image_index = static_cast<int>(images.size());
        auto format = srgb_indices.contains(image_index) ? GL_SRGB8_ALPHA8 : GL_RGBA8;
        glTextureStorage2D(texture.id(), GetMipmapLevels(width, height), format, width, height);
        glTextureSubImage2D(texture.id(), 0, 0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, pbuffer->data());
        glGenerateTextureMipmap(texture.id());

        images.push_back(std::move(texture));
    }

    return images;
}

std::vector<int> LoadTextures(const tinygltf::Model& input) {
    // TODO: Sampler

    std::vector<int> textures;
    textures.reserve(input.textures.size());
    for (const auto& texture : input.textures) {
        textures.push_back(texture.source);
    }
    return textures;
}

std::vector<Node> LoadNodes(const tinygltf::Model& input) {
    std::vector<Node> nodes;
    nodes.reserve(input.nodes.size());
    for (const auto& node : input.nodes) {
        Node n;
        n.children = node.children;
        n.mesh = node.mesh;

        if (node.translation.size() == 3) {
            n.translation = glm::make_vec3(node.translation.data());
        }
        if (node.rotation.size() == 4) {
            n.rotation = glm::make_quat(node.rotation.data());
#ifndef GLM_FORCE_QUAT_DATA_XYZW
            // https://github.com/g-truc/glm/commit/820a2c0e625f26000c688d841836bb10483be34d
            n.rotation = { n.rotation[3], n.rotation[0], n.rotation[1], n.rotation[2] };
#endif
        }
        if (node.scale.size() == 3) {
            n.scale = glm::make_vec3(node.scale.data());
        }
        if (node.matrix.size() == 16) {
            n.matrix = glm::make_mat4x4(node.matrix.data());
        }
        nodes.push_back(n);
    }
    return nodes;
}

struct MeshesData{
    std::vector<GltfModel::Mesh> meshes;
    std::vector<Vertex> vertex_buffer;
    std::vector<int> index_buffer;
} LoadMeshes(const tinygltf::Model& input) {
    MeshesData res;
    res.meshes.reserve(input.meshes.size());

    for (const auto& mesh : input.meshes) {
        GltfModel::Mesh m;
        m.primitives.reserve(mesh.primitives.size());
        for (const auto& primitive : mesh.primitives) {
            auto base_index = static_cast<uint32_t>(res.index_buffer.size());
            auto vertex_start = static_cast<uint32_t>(res.vertex_buffer.size());

            // Vertices
            const float* position_buffer = nullptr;
            const float* normal_buffer = nullptr;
            const float* uv_buffer = nullptr;
            const float* tangent_buffer = nullptr;
            size_t vertex_count = 0;

            for (auto [pp_buffer, key] : { 
                    std::make_pair(&position_buffer, "POSITION"),
                    std::make_pair(&normal_buffer, "NORMAL"),
                    std::make_pair(&uv_buffer, "TEXCOORD_0"),
                    std::make_pair(&tangent_buffer, "TANGENT"), }) {
                if (auto itr = primitive.attributes.find(key); itr != primitive.attributes.end()) {
                    const auto& accessor = input.accessors[itr->second];
                    ASSERT(accessor.componentType == TINYGLTF_PARAMETER_TYPE_FLOAT);
                    const auto& view = input.bufferViews[accessor.bufferView];
                    *pp_buffer = reinterpret_cast<const float*>(&(input.buffers[view.buffer].data[accessor.byteOffset + view.byteOffset]));
                    vertex_count = accessor.count;
                }
            }

            ASSERT(position_buffer);
            ASSERT(normal_buffer);

            for (size_t v = 0; v < vertex_count; v++) {
                Vertex vert{};
                vert.position = glm::make_vec3(&position_buffer[v * 3]);
                vert.normal = glm::make_vec3(&normal_buffer[v * 3]);
                vert.uv = uv_buffer ? glm::make_vec2(&uv_buffer[v * 2]) : glm::vec2(0.0f);
                vert.uv.y = 1.0f - vert.uv.y; // flip y
                if (tangent_buffer) {
                    // TODO: flip y
                    vert.tangent = glm::normalize(glm::make_vec3(&tangent_buffer[v * 3]));
                }
                else {
                    const glm::vec3 up(0.0f, 1.0f, 0.0f);
                    const glm::vec3 right(1.0f, 0.0f, 0.0f);
                    // TODO: calculate tangents
                    vert.tangent = glm::normalize(glm::cross(vert.normal, 
                        abs(glm::dot(vert.normal, up)) > 0.01f ? right : up));
                }
                res.vertex_buffer.push_back(vert);
            }

            // Indices
            const auto& accessor = input.accessors[primitive.indices];
            const auto& view = input.bufferViews[accessor.bufferView];
            const auto& buffer = input.buffers[view.buffer];

            auto copy_data = [&accessor, &res, vertex_start](auto buf) {
                for (size_t index = 0; index < accessor.count; index++) {
                    res.index_buffer.push_back(buf[index] + vertex_start);
                }
            };
            auto p_data = &buffer.data[accessor.byteOffset + view.byteOffset];
            switch (accessor.componentType) {
            case TINYGLTF_PARAMETER_TYPE_UNSIGNED_INT:
                copy_data(reinterpret_cast<const uint32_t*>(p_data));
                break;
            case TINYGLTF_PARAMETER_TYPE_UNSIGNED_SHORT:
                copy_data(reinterpret_cast<const uint16_t*>(p_data));
                break;
            case TINYGLTF_PARAMETER_TYPE_UNSIGNED_BYTE:
                copy_data(reinterpret_cast<const uint8_t*>(p_data));
                break;
            default:
                ASSERT(false);
            }

            Primitive p;
            p.material_index = primitive.material;
            p.mode = primitive.mode;
            p.first_index = base_index;
            p.count = static_cast<uint32_t>(accessor.count);
            m.primitives.push_back(p);
        }
        res.meshes.push_back(std::move(m));
    }

    return res;
}

GltfModel::GltfModel(const std::filesystem::path& path) {
    tinygltf::TinyGLTF gltf_ctx;
    tinygltf::Model model;
    std::string err;
    std::string warn;
    std::string ext = path.extension().string();

    LOG_INFO("Loading {}", path.string());

    bool ret = false;
    if (ext.compare(".glb") == 0) {
        ret = gltf_ctx.LoadBinaryFromFile(&model, &err, &warn, path.string());
    }
    else {
        ret = gltf_ctx.LoadASCIIFromFile(&model, &err, &warn, path.string());
    }

    if (!warn.empty()) {
        LOG_WARN("{}", warn);
    }

    if (!err.empty()) {
        LOG_ERROR("{}", err);
    }

    ASSERT(ret);

    materials_ = LoadMaterials(model);
    textures_ = LoadTextures(model);

    std::unordered_set<int> srgb_indices;
    for (const auto& m : materials_) {
        for (auto tex_i : { m.base_color_index , m.emissive_index }) {
            if (tex_i < 0)
                continue;
            auto img_i = textures_[tex_i];
            if (img_i < 0)
                continue;
            srgb_indices.insert(img_i);
        }
    }

    images_ = LoadImages(model, srgb_indices);
    nodes_ = LoadNodes(model);
    ASSERT(model.scenes.size() > 0);
    scene_nodes_ = model.scenes[0].nodes;
    auto meshes_data = LoadMeshes(model);
    const auto& vertex_buffer = meshes_data.vertex_buffer;
    const auto& index_buffer = meshes_data.index_buffer;

    vao_.Create();
    vbo_.Create();
    ebo_.Create();
    constexpr int kStride = 11;
    glNamedBufferStorage(vbo_.id(), sizeof(vertex_buffer[0]) * vertex_buffer.size(), vertex_buffer.data(), 0);
    glNamedBufferStorage(ebo_.id(), sizeof(index_buffer[0]) * index_buffer.size(), index_buffer.data(), 0);
    glBindVertexArray(vao_.id());
    glBindBuffer(GL_ARRAY_BUFFER, vbo_.id());
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo_.id());
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, kStride * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, kStride * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, kStride * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, kStride * sizeof(float), (void*)(8 * sizeof(float)));
    glEnableVertexAttribArray(3);

    meshes_ = std::move(meshes_data.meshes);
}

void GltfModel::UpdateMatrices(const glm::mat4& m) {
    matrices_.resize(nodes_.size());
    for (auto node : scene_nodes_) {
        UpdateMatrices(m, node);
    }
}

void GltfModel::UpdateMatrices(const glm::mat4& m, int i) {
    matrices_[i] = m * nodes_[i].LocalMatrix();
    for (auto child : nodes_[i].children) {
        UpdateMatrices(matrices_[i], child);
    }
}

void DrawPrimitive(const Primitive& primitive) {
    glDrawElements(primitive.mode, primitive.count, GL_UNSIGNED_INT,
        reinterpret_cast<void*>(sizeof(uint32_t) * primitive.first_index));
}

void GltfModel::RenderToGBuffer(const glm::mat4& view_projection) const {
    glBindVertexArray(vao_.id());
    for (size_t i = 0; i < nodes_.size(); ++i) {
        auto mesh_i = nodes_[i].mesh;
        if (mesh_i < 0)
            continue;
        for (const auto& primitive : meshes_[mesh_i].primitives) {
            ::Material mesh_object_material;
            const auto& material = materials_[primitive.material_index];
            mesh_object_material.albedo_factor = material.base_color_factor;
            mesh_object_material.albedo_texture = GetTexture(material.base_color_index);
            mesh_object_material.metallic_factor = 1.0f;
            mesh_object_material.roughness_factor = 1.0f;
            mesh_object_material.ao_factor = material.ao_factor;
            mesh_object_material.normal_texture = GetTexture(material.normal_index);
            mesh_object_material.metallic_roughness_texture = GetTexture(material.metallic_roughness_index);
            mesh_object_material.ao_texture = GetTexture(material.ao_index);

            GBufferRenderer::Instance().Setup(matrices_[i], view_projection, mesh_object_material);
            DrawPrimitive(primitive);
        }
    }
}

void GltfModel::RenderToShadowMap(const glm::mat4& light_view_projection) const {
    glBindVertexArray(vao_.id());
    for (size_t i = 0; i < nodes_.size(); ++i) {
        auto mesh_i = nodes_[i].mesh;
        if (mesh_i < 0)
            continue;
        ShadowMapRenderer::Instance().Setup(matrices_[i], light_view_projection);
        for (const auto& primitive : meshes_[mesh_i].primitives) {
            DrawPrimitive(primitive);
        }
    }
}
