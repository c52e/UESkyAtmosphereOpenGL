#include "SMAA.h"

#include <sstream>

#include "ScreenRectangle.h"
#include "Samplers.h"

#include "SMAA/AreaTex.h"
#include "SMAA/SearchTex.h"
#include "log.h"

SMAA::SMAA(int width, int height, SMAAOption option)
	: option_(option) {
	auto load = [width, height, option](const char* path) {
		return [path, width, height, option]() {
			auto src = ReadWithPreprocessor(path);
			std::stringstream ss;
			ss << "#version 460\n" << "#define SMAA_RT_METRICS vec4("
				<< std::to_string(1.0 / width) << ","
				<< std::to_string(1.0 / height) << ","
				<< std::to_string(width) << ","
				<< std::to_string(height) << ")\n";
			ss << "#define " << magic_enum::enum_name(option) << "\n";
			auto common = ss.str();
			auto vertex = common + "#define VERTEX\n" + src;
			auto fragment = common + "#define FRAGMENT\n" + src;
			return GLProgram(vertex.c_str(), fragment.c_str());
		};
	};;

	edge_detection_ = load("../shaders/Base/SMAA/EdgeDetection.glsl");
	blending_weight_calculation_ = load("../shaders/Base/SMAA/BlendingWeightCalculation.glsl");
	neighborhood_blending_ = load("../shaders/Base/SMAA/NeighborhoodBlending.glsl");

	auto create_tex = [](const unsigned char data[], int width, int height, int element_bytes) {
		ASSERT(element_bytes == 1 || element_bytes == 2);

		GLTexture res;
		res.Create(GL_TEXTURE_2D);

		GLenum internalformat[]{0, GL_R8, GL_RG8};
		GLenum formats[]{0, GL_RED, GL_RG};
		glTextureStorage2D(res.id(), 1, internalformat[element_bytes], width, height);

		std::vector<unsigned char> processed(data, data + width * height * element_bytes);
		glTextureSubImage2D(res.id(), 0, 0, 0, width, height, formats[element_bytes], GL_UNSIGNED_BYTE, processed.data());

		return res;
	};
	area_tex_ = create_tex(areaTexBytes, AREATEX_WIDTH, AREATEX_HEIGHT, 2);
	search_tex_ = create_tex(searchTexBytes, SEARCHTEX_WIDTH, SEARCHTEX_HEIGHT, 1);


	framebuffer_.Create();
	edges_tex_.Create(GL_TEXTURE_2D);
	blend_tex_.Create(GL_TEXTURE_2D);
	output_tex_.Create(GL_TEXTURE_2D);
	stencil_tex_.Create(GL_TEXTURE_2D);
	glTextureStorage2D(edges_tex_.id(), 1, GL_RGBA8, width, height);
	glTextureStorage2D(blend_tex_.id(), 1, GL_RGBA8, width, height);
	glTextureStorage2D(output_tex_.id(), 1, GL_SRGB8_ALPHA8, width, height);
	glTextureStorage2D(stencil_tex_.id(), 1, GL_STENCIL_INDEX8, width, height);

	glNamedFramebufferTexture(framebuffer_.id(), GL_COLOR_ATTACHMENT0, edges_tex_.id(), 0);
	GLenum attachments[]{ GL_COLOR_ATTACHMENT0 };
	glNamedFramebufferDrawBuffers(framebuffer_.id(), GLsizei(std::size(attachments)), attachments);
	glNamedFramebufferTexture(framebuffer_.id(), GL_STENCIL_ATTACHMENT, stencil_tex_.id(), 0);
	if (glCheckNamedFramebufferStatus(framebuffer_.id(), GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
		throw std::runtime_error("Framebuffer Incomplete");
	}

}

void SMAA::DoSMAA(GLuint input) {
	PERF_MARKER("SMAA")
	glBindFramebuffer(GL_FRAMEBUFFER, framebuffer_.id());
	if (option_ == SMAAOption::OFF) {
		glNamedFramebufferTexture(framebuffer_.id(), GL_COLOR_ATTACHMENT0, output_tex_.id(), 0);
		TextureVisualizer::Instance().VisualizeTexture(input);
		return;
	}

	EdgesDetectionPass(input);
	BlendingWeightsCalculationPass();
	NeighborhoodBlendingPass(input);
}

void SMAA::EdgesDetectionPass(GLuint input) {
	PERF_MARKER("EdgesDetectionPass")
	glNamedFramebufferTexture(framebuffer_.id(), GL_COLOR_ATTACHMENT0, edges_tex_.id(), 0);

	const float kBlack[] = { 0.f, 0.f, 0.f, 0.f };
	glClearBufferfv(GL_COLOR, 0, kBlack);
	GLint kZero = 0;
	glClearBufferiv(GL_STENCIL, 0, &kZero);

	glEnable(GL_STENCIL_TEST);
	glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
	glStencilMask(0xff);
	glStencilFunc(GL_ALWAYS, 1, 0xff);

	GLBindTextures({ input });
	GLBindSamplers({ Samplers::GetLinearNoMipmapClampToEdge() });
	glUseProgram(edge_detection_.id());
	ScreenRectangle::Instance().Draw();
}

void SMAA::BlendingWeightsCalculationPass() {
	PERF_MARKER("BlendingWeightsCalculationPass")
	glNamedFramebufferTexture(framebuffer_.id(), GL_COLOR_ATTACHMENT0, blend_tex_.id(), 0);
	const float kBlack[] = { 0.f, 0.f, 0.f, 0.f };
	glClearBufferfv(GL_COLOR, 0, kBlack);

	glStencilFunc(GL_EQUAL, 1, 0xff);
	glStencilMask(0x00);

	GLBindTextures({ edges_tex_.id(),
		area_tex_.id(),
		search_tex_.id() });
	GLBindSamplers({ Samplers::GetLinearNoMipmapClampToEdge(),
		Samplers::GetLinearNoMipmapClampToEdge(),
		Samplers::GetLinearNoMipmapClampToEdge() });
	glUseProgram(blending_weight_calculation_.id());
	ScreenRectangle::Instance().Draw();

	glStencilMask(0xff);
	glDisable(GL_STENCIL_TEST);
}

void SMAA::NeighborhoodBlendingPass(GLuint input) {
	PERF_MARKER("NeighborhoodBlendingPass")
	glNamedFramebufferTexture(framebuffer_.id(), GL_COLOR_ATTACHMENT0, output_tex_.id(), 0);
	GLBindTextures({ input,
		blend_tex_.id() });
	GLBindSamplers({ Samplers::GetLinearNoMipmapClampToEdge(),
		Samplers::GetLinearNoMipmapClampToEdge() });
	glUseProgram(neighborhood_blending_.id());
	ScreenRectangle::Instance().Draw();
}
