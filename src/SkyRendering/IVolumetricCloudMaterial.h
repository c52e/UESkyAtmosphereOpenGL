#pragma once

#include <glm/glm.hpp>

#include "gl.hpp"
#include "Camera.h"
#include "Serialization.h"

class IVolumetricCloudMaterial : public ISerializable {
public:
    static constexpr GLuint kMaterialTextureUnitBegin = 7;

    virtual ~IVolumetricCloudMaterial() = default;

    virtual std::string ShaderSrc() = 0;

    struct UpdateParam {
        glm::vec2 viewport;
        const Camera& camera;
        glm::dvec2 offset_from_first;
        glm::dmat4 model;
    };

    virtual void Update(const UpdateParam& param, glm::vec2& additional_delta) {}

    virtual void Bind() = 0;

    virtual void DrawGUI() = 0;

    virtual float GetSigmaTMax() = 0;

    FIELD_DECLARATION_BEGIN(ISerializable)
    FIELD_DECLARATION_END()
};

HAS_SUBCLASS(IVolumetricCloudMaterial)
