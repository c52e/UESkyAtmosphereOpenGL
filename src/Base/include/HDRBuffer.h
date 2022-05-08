#pragma once

#include <glm/glm.hpp>
#include <magic_enum.hpp>

#include "gl.hpp"
#include "Singleton.h"
#include "GLReloadableProgram.h"
#include "Serialization.h"

enum class ToneMapping {
	CEToneMapping,
	ACESToneMapping
};

struct PostProcessParameters : public ISerializable {
	ToneMapping tone_mapping = ToneMapping::ACESToneMapping;
	float bloom_min_luminance = 1.0f;
	float bloom_max_delta_luminance = 20.0f;
	float bloom_filter_width = 0.03f;
	float bloom_intensity = 0.1f;
	float exposure = 10.f;
	bool dither_color_enable = true;

	FIELD_DECLARATION_BEGIN(ISerializable)
		FIELD_DECLARE(tone_mapping)
		FIELD_DECLARE(bloom_min_luminance)
		FIELD_DECLARE(bloom_max_delta_luminance)
		FIELD_DECLARE(bloom_filter_width)
		FIELD_DECLARE(bloom_intensity)
		FIELD_DECLARE(exposure)
		FIELD_DECLARE(dither_color_enable)
	FIELD_DECLARATION_END()
};

class HDRBuffer {
public:
	HDRBuffer(int width, int height);

	GLuint id() const { return framebuffers_[0]; }
	GLuint texture_id() const { return textures_[0]; }
	void DoPostProcess(GLuint target_framebuffer_id, const PostProcessParameters& params);

private:
	GLFramebuffers<3> framebuffers_;
	GLTextures<3> textures_;

	class PostProcessRenderer : public Singleton<PostProcessRenderer> {
	public:
		friend Singleton<PostProcessRenderer>;

		void Render(GLuint target_framebuffer_id, const HDRBuffer& hdrbuffer, const PostProcessParameters& params);
	private:
		PostProcessRenderer();
		GLReloadableProgram extract_;
		GLReloadableProgram pass1_;
		GLReloadableProgram pass2_[magic_enum::enum_count<ToneMapping>()][2];
	};
};

