#include "Mesh.h"

#include <cmath>
#include <fstream>
#include <filesystem>

#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <tiny_gltf.h>

#include "log.h"

inline glm::vec3 GetAnyOrthogonalVector(const glm::vec3& N) {
    if (glm::abs(N.z) < 1e-6f)
        return { 0,0,1 };
    return glm::normalize(glm::vec3(1, 1, (-N.x - N.y) / N.z));
}

Mesh::Mesh(const MeshVertices& vertices) {
    ASSERT(vertices.mode == GL_TRIANGLES
        || vertices.mode == GL_TRIANGLE_STRIP
        || vertices.mode == GL_TRIANGLE_FAN);
    ASSERT(vertices.positions.size() == vertices.normals.size());

    auto vertices_num = vertices.positions.size();
    mode = vertices.mode;
    count = static_cast<GLsizei>(vertices.indices.size());
    constexpr int stride = 11;

    vao_.Create();
    vbo_.Create();
    ebo_.Create();
    std::vector<float> data;
    data.reserve(vertices_num * stride);
    for (int i = 0; i < vertices_num; ++i) {
        data.push_back(vertices.positions[i][0]);
        data.push_back(vertices.positions[i][1]);
        data.push_back(vertices.positions[i][2]);
        data.push_back(vertices.normals[i][0]);
        data.push_back(vertices.normals[i][1]);
        data.push_back(vertices.normals[i][2]);
        if (!vertices.uvs.empty()) {
            data.push_back(vertices.uvs[i][0]);
            data.push_back(vertices.uvs[i][1]);
            if (vertices.tangents.empty()) {
               // TODO: calculate tangents
                throw std::runtime_error("Tangents is empty");
            }
            data.push_back(vertices.tangents[i][0]);
            data.push_back(vertices.tangents[i][1]);
            data.push_back(vertices.tangents[i][2]);
        }
        else {
            data.push_back(0);
            data.push_back(0);
            auto tangent = GetAnyOrthogonalVector(glm::make_vec3(vertices.normals[i].data()));
            data.push_back(tangent[0]);
            data.push_back(tangent[1]);
            data.push_back(tangent[2]);
        }
    }
    glNamedBufferStorage(vbo_.id(), sizeof(data[0]) * data.size(), data.data(), 0);
    glNamedBufferStorage(ebo_.id(), sizeof(vertices.indices[0]) * vertices.indices.size(),
        vertices.indices.data(), 0);
    glBindVertexArray(vao_.id());
    glBindBuffer(GL_ARRAY_BUFFER, vbo_.id());
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo_.id());
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, stride * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, stride * sizeof(float), (void*)(8 * sizeof(float)));
    glEnableVertexAttribArray(3);
}

void Mesh::Draw() const {
    glBindVertexArray(vao_.id());
    glDrawElements(mode, count, GL_UNSIGNED_INT, 0);
}

MeshVertices CreatePyramid() {
    MeshVertices vertices;
    vertices.mode = GL_TRIANGLES;
    vertices.positions = {
        {0.f, 1.f, 0.f}, {1.f, 0.f, 1.f}, {1.f, 0.f, -1.f},
        {0.f, 1.f, 0.f}, {1.f, 0.f, -1.f}, {-1.f, 0.f, -1.f},
        {0.f, 1.f, 0.f}, {-1.f, 0.f, -1.f}, {-1.f, 0.f, 1.f},
        {0.f, 1.f, 0.f}, {-1.f, 0.f, 1.f}, {1.f, 0.f, 1.f}};
    vertices.normals = {
        {1.f, 1.f, 0.f}, {1.f, 1.f, 0.f}, {1.f, 1.f, 0.f},
        {0.f, 1.f, -1.f}, {0.f, 1.f, -1.f}, {0.f, 1.f, -1.f},
        {-1.f, 1.f, 0.f}, {-1.f, 1.f, 0.f}, {-1.f, 1.f, 0.f},
        {0.f, 1.f, 1.f}, {0.f, 1.f, 1.f}, {0.f, 1.f, 1.f}};
    vertices.indices = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
    return vertices;
}

MeshVertices CreateSphere() {
    // https://learnopengl.com/code_viewer_gh.php?code=src/6.pbr/1.2.lighting_textured/lighting_textured.cpp
    MeshVertices vertices;
    vertices.mode = GL_TRIANGLE_STRIP;

    auto pi = static_cast<float>(PI);
    const unsigned int X_SEGMENTS = 64;
    const unsigned int Y_SEGMENTS = 64;
    for (unsigned int y = 0; y <= Y_SEGMENTS; ++y) {
        for (unsigned int x = 0; x <= X_SEGMENTS; ++x) {
            float xSegment = (float)x / (float)X_SEGMENTS;
            float ySegment = (float)y / (float)Y_SEGMENTS;
            float xPos = cos(xSegment * 2.0f * pi) * sin(ySegment * pi);
            float yPos = cos(ySegment * pi);
            float zPos = -sin(xSegment * 2.0f * pi) * sin(ySegment * pi);

            float xTan = -sin(xSegment * 2.0f * pi);
            float yTan = 0;
            float zTan = -cos(xSegment * 2.0f * pi);

            vertices.positions.push_back({ xPos, yPos, zPos });
            vertices.normals.push_back({ xPos, yPos, zPos });
            vertices.tangents.push_back({ xTan , yTan, zTan });
            vertices.uvs.push_back({ xSegment, 1.0f - ySegment });
        }
    }

    bool oddRow = false;
    for (int y = 0; y < Y_SEGMENTS; ++y) {
        if (!oddRow) {
            for (int x = 0; x <= X_SEGMENTS; ++x) {
                vertices.indices.push_back(y * (X_SEGMENTS + 1) + x);
                vertices.indices.push_back((y + 1) * (X_SEGMENTS + 1) + x);
            }
        }
        else {
            for (int x = X_SEGMENTS; x >= 0; --x) {
                vertices.indices.push_back((y + 1) * (X_SEGMENTS + 1) + x);
                vertices.indices.push_back(y * (X_SEGMENTS + 1) + x);
            }
        }
        oddRow = !oddRow;
    }
    return vertices;
}

namespace fs = std::filesystem;

bool LoadModel(tinygltf::Model& model, const char* path) {
    tinygltf::TinyGLTF gltf_ctx;
    std::string err;
    std::string warn;
    std::string ext = fs::path(path).extension().string();

    LOG_INFO("Loading {}", path);

    bool ret = false;
    if (ext.compare(".glb") == 0) {
        ret = gltf_ctx.LoadBinaryFromFile(&model, &err, &warn, path);
    }
    else {
        ret = gltf_ctx.LoadASCIIFromFile(&model, &err, &warn, path);
    }

    if (!warn.empty()) {
        LOG_WARN("{}", warn);
    }

    if (!err.empty()) {
        LOG_ERROR("{}", err);
    }

    return ret;
}

MeshVertices ReadMeshFile(const char* path) {
    tinygltf::Model m;
    if (!LoadModel(m, path)) {
        throw std::runtime_error("Failed to parse glTF\n");
    }

    MeshVertices vertices;
    vertices.mode = GL_TRIANGLES;

    auto node_i = m.scenes[0].nodes[0];
    auto mesh_i = m.nodes[node_i].mesh;
    auto& mesh = m.meshes[mesh_i];
    auto& primitive = mesh.primitives[0];

    ASSERT(primitive.indices == 0);
    ASSERT(primitive.attributes["POSITION"] == 1);
    ASSERT(primitive.attributes["NORMAL"] == 2);

    ASSERT(m.accessors[0].componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT);
    ASSERT(m.accessors[0].type == TINYGLTF_TYPE_SCALAR);
    ASSERT(m.accessors[1].componentType == TINYGLTF_COMPONENT_TYPE_FLOAT);
    ASSERT(m.accessors[1].type == TINYGLTF_TYPE_VEC3);
    ASSERT(m.accessors[2].componentType == TINYGLTF_COMPONENT_TYPE_FLOAT);
    ASSERT(m.accessors[2].type == TINYGLTF_TYPE_VEC3);

    auto copydata = [&m](int i, auto& v) {
        auto& accessor = m.accessors[i];
        v.resize(accessor.count);
        auto offset = accessor.byteOffset;
        auto size = sizeof(v[0]) * v.size();

        auto& bufferView = m.bufferViews[accessor.bufferView];
        auto& buffer = m.buffers[bufferView.buffer];
        memcpy(v.data(), buffer.data.data() + bufferView.byteOffset + offset, size);
    };
    copydata(0, vertices.indices);
    copydata(1, vertices.positions);
    copydata(2, vertices.normals);

    return vertices;
}
