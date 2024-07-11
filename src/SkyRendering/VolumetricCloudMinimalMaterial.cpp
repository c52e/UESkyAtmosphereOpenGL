#include "VolumetricCloudMinimalMaterial.h"

#include <imgui.h>

#include "GLProgram.h"

struct VolumetricCloudMinimalMaterial::BufferData {
	glm::vec3 padding;
	float uDensity;
};

VolumetricCloudMinimalMaterial::VolumetricCloudMinimalMaterial() {
	buffer_.Create();
	glNamedBufferStorage(buffer_.id(), sizeof(BufferData), NULL, GL_DYNAMIC_STORAGE_BIT);
}

std::string VolumetricCloudMinimalMaterial::ShaderSrc() {
	return ReadWithPreprocessor("../shaders/SkyRendering/VolumetricCloudMaterialMinimal.glsl");
}

void VolumetricCloudMinimalMaterial::Bind() {
	BufferData buffer;
	buffer.uDensity = density_;

	glNamedBufferSubData(buffer_.id(), 0, sizeof(buffer), &buffer);
	glBindBufferBase(GL_UNIFORM_BUFFER, 3, buffer_.id());
}

float VolumetricCloudMinimalMaterial::GetSigmaTMax() {
	return density_;
}

void VolumetricCloudMinimalMaterial::DrawGUI() {
	ImGui::SliderFloat("Density", &density_, 0.0f, 1.0f);
}
