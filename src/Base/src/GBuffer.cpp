#include "GBuffer.h"

#include <stdexcept>
#include <array>

#include <glm/gtc/type_ptr.hpp>

#include "Textures.h"
#include "Samplers.h"
#include "Utils.h"

GBuffer::GBuffer(int width, int height) 
	: dim_{ width, height } {
	framebuffer_.Create();
	albedo_.Create(GL_TEXTURE_2D);
	normal_.Create(GL_TEXTURE_2D);
	orm_.Create(GL_TEXTURE_2D);
	depth_.Create(GL_TEXTURE_2D);

	glTextureStorage2D(albedo_.id(), 1, GL_SRGB8_ALPHA8, width, height);
	glTextureStorage2D(normal_.id(), 1, GL_RGBA16_SNORM, width, height);
	glTextureStorage2D(orm_.id(), 1, GL_RGBA16, width, height);
	glTextureStorage2D(depth_.id(), 1, GL_DEPTH_COMPONENT32F, width, height);

	glNamedFramebufferTexture(framebuffer_.id(), GL_COLOR_ATTACHMENT0, albedo_.id(), 0);
	glNamedFramebufferTexture(framebuffer_.id(), GL_COLOR_ATTACHMENT1, normal_.id(), 0);
	glNamedFramebufferTexture(framebuffer_.id(), GL_COLOR_ATTACHMENT2, orm_.id(), 0);
	glNamedFramebufferTexture(framebuffer_.id(), GL_DEPTH_ATTACHMENT, depth_.id(), 0);

	GLenum attachments[]{ GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2 };
	constexpr GLsizei attachments_num = sizeof(attachments) / sizeof(attachments[0]);
	glNamedFramebufferDrawBuffers(framebuffer_.id(), attachments_num, attachments);

	if (glCheckNamedFramebufferStatus(framebuffer_.id(), GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
		throw std::runtime_error("GBuffer Incomplete");
	}

	pixel_visibility_.Create(GL_TEXTURE_2D);
	glTextureStorage2D(pixel_visibility_.id(), 1, kVisiblityTextureFormat, width, height);
}

GBufferRenderer::GBufferRenderer() {
	program_ = []() {
		auto src = ReadFile("../shaders/Base/GBuffer.glsl");
		auto vertex = Replace(src, "VERTEX", "1");
		auto fragment = Replace(src, "FRAGMENT", "1");
		return GLProgram(vertex.c_str(), fragment.c_str());
	};
}

void GBufferRenderer::Setup(const glm::mat4 model_matrix, const glm::mat4& view_projection_matrix, const Material& material) {
	glUseProgram(program_.id());
	auto mvp = view_projection_matrix * model_matrix;
	glUniformMatrix4fv(0, 1, GL_FALSE, glm::value_ptr(mvp));
	auto normal_matrix = glm::mat3(glm::transpose(glm::inverse(model_matrix)));
	glUniformMatrix3fv(1, 1, GL_FALSE, glm::value_ptr(normal_matrix));
	glUniform3fv(2, 1, glm::value_ptr(material.albedo_factor));
	glUniform2f(3, material.metallic_factor, material.roughness_factor);
	glUniform1f(4, material.ao_factor);
	GLBindTextures({ Textures::Instance().white(material.albedo_texture),
					Textures::Instance().normal(material.normal_texture),
					Textures::Instance().white(material.metallic_roughness_texture),
					Textures::Instance().white(material.ao_texture) });
	GLBindSamplers({ Samplers::GetAnisotropySampler(Samplers::Wrap::REPEAT),
					Samplers::GetAnisotropySampler(Samplers::Wrap::REPEAT),
					Samplers::GetAnisotropySampler(Samplers::Wrap::REPEAT),
					Samplers::GetAnisotropySampler(Samplers::Wrap::REPEAT) });
}
