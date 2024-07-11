#pragma once

#include <glm/glm.hpp>

#include "gl.hpp"
#include "Singleton.h"
#include "GLReloadableProgram.h"
#include "MeshObject.h"

class GBuffer {
public:
	static constexpr GLenum kVisiblityTextureFormat = GL_R16;

	GBuffer(int width, int height);

	glm::ivec2 dim() const { return dim_; }

	GLuint id() const { return framebuffer_.id(); }
	GLuint albedo() const { return albedo_.id(); }
	GLuint normal() const { return normal_.id(); }
	GLuint orm() const { return orm_.id(); }
	GLuint depth() const { return depth_.id(); }

	GLuint pixel_visibility() const { return pixel_visibility_.id(); }

private:
	glm::ivec2 dim_;

	GLFramebuffer framebuffer_;
	GLTexture albedo_;
	GLTexture normal_;
	GLTexture orm_;
	GLTexture depth_;

	GLTexture pixel_visibility_;
};

inline void Clear(const GBuffer& gbuffer) {
	static const float black[] = { 0.f, 0.f, 0.f, 0.f };
	glClearNamedFramebufferfv(gbuffer.id(), GL_COLOR, 0, black);
	glClearNamedFramebufferfv(gbuffer.id(), GL_COLOR, 1, black);
	glClearNamedFramebufferfv(gbuffer.id(), GL_COLOR, 2, black);
	
	glClearNamedFramebufferfv(gbuffer.id(), GL_DEPTH, 0, black);

	float one = 1.0f;
	glClearTexImage(gbuffer.pixel_visibility(), 0, GL_RED, GL_FLOAT, &one);
}


class GBufferRenderer : public Singleton<GBufferRenderer> {
public:
	friend Singleton<GBufferRenderer>;

	void Setup(const glm::mat4 model_matrix, const glm::mat4& view_projection_matrix, const Material& material);

private:
	GBufferRenderer();

	GLReloadableProgram program_;
};
