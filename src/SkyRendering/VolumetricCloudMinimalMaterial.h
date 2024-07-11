#pragma once

#include "IVolumetricCloudMaterial.h"

class VolumetricCloudMinimalMaterial : public IVolumetricCloudMaterial {
public:
    VolumetricCloudMinimalMaterial();

    std::string ShaderSrc() override;

    void Bind() override;

    float GetSigmaTMax() override;

    void DrawGUI() override;

private:
    FIELD_DECLARATION_BEGIN(ISerializable)
        FIELD_DECLARE(density_)
    FIELD_DECLARATION_END()

    struct BufferData;

    GLBuffer buffer_;

    float density_ = 0.5f;
};
