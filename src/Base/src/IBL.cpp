#include "IBL.h"

#include "Samplers.h"

IBL::IBL() {
    env_radiance_sh_program_ = []() {
        auto src = ReadWithPreprocessor("../shaders/Base/EnvRadianceSH.comp");
        return GLProgram(src.c_str());
    };

    env_radiance_sh_buffer_.Create();
    glNamedBufferStorage(env_radiance_sh_buffer_.id(), sizeof(glm::vec4) * 9, NULL, GL_DYNAMIC_STORAGE_BIT);

    prefilter_radiance_program_ = {
        "../shaders/Base/PrefilterRadiance.comp",
        {{8, 4, 1}, {8, 8, 1}, {16, 4, 1}, {16, 8, 1}}
    };

    prefiltered_radiance_.Create(GL_TEXTURE_CUBE_MAP);
    glTextureStorage2D(prefiltered_radiance_.id(), kRoughnessCount, kPrefilteredRadianceFormat, kPrefilteredRadianceResolution, kPrefilteredRadianceResolution);
}

void IBL::Precompute(GLuint environment_radiance_texture) {
    GLBindTextures({ environment_radiance_texture });
    GLBindSamplers({ Samplers::GetAnisotropySampler(Samplers::Wrap::CLAMP_TO_EDGE) });
    {
        PERF_MARKER("Environment Radiance SH")
        glUseProgram(env_radiance_sh_program_.id());
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, env_radiance_sh_buffer_.id());
        glDispatchCompute(9, 1, 1);
    }
    {
        PERF_MARKER("Prefilter Radiance")
        glUseProgram(prefilter_radiance_program_.id());
        for (int i = 0, w = kPrefilteredRadianceResolution; i < kRoughnessCount; ++i, w >>= 1) {
            glBindImageTexture(0, prefiltered_radiance_.id(), i, GL_TRUE, 0, GL_WRITE_ONLY, kPrefilteredRadianceFormat);
            glUniform1f(0, float(i) / float(kRoughnessCount - 1));
            prefilter_radiance_program_.Dispatch({ w, w, 6 });
        }
    }

    glMemoryBarrier(GL_UNIFORM_BARRIER_BIT | GL_TEXTURE_FETCH_BARRIER_BIT);
}
