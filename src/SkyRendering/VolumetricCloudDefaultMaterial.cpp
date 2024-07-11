#include "VolumetricCloudDefaultMaterial.h"

#include <imgui.h>

#include "ImageLoader.h"
#include "Samplers.h"
#include "ImGuiExt.h"

struct SampleInfo {
	glm::vec2 bias;
	float frequency;
	float k_lod;
};

struct VolumetricCloudDefaultMaterialCommonBufferData {
	SampleInfo uCloudMapSampleInfo;
	SampleInfo uDetailSampleInfo;
	SampleInfo uDisplacementSampleInfo;
	glm::vec3 padding0;
	float uDensity;
};

float CalKLod(const TextureWithInfo& tex, glm::vec2 viewport, const Camera& camera) {
	auto max_width = static_cast<float>(glm::max(tex.x, glm::max(tex.y, tex.z)));
	return max_width * camera.tangent_half_fovy() * 2.0f / (tex.repeat_size * static_cast<float>(viewport.y));
}

VolumetricCloudDefaultMaterialCommon::VolumetricCloudDefaultMaterialCommon(int cloud_map_size, int detail_size, int displacement_size) {
	{
		const int w = cloud_map_size;
		cloud_map_.texture.x = cloud_map_.texture.y = w;
		cloud_map_.texture.tex.Create(GL_TEXTURE_2D);
		glTextureStorage2D(cloud_map_.texture.id(), GetMipmapLevels(w, w), GL_RG8, w, w);
		cloud_map_.texture.repeat_size = 18.99f;

		cloud_map_.program = {
			"../shaders/SkyRendering/NoiseGen.comp",
			{{16, 8}, {32, 16}, {32, 32}, {8, 8}, {16, 16}, {8, 4}},
			[](const std::string& src) { return Replace(src, "CLOUD_MAP_GEN", "1"); },
			"CLOUD_MAP"
		};
	}
	{
		const int w = detail_size;
		detail_.texture.x = detail_.texture.y = detail_.texture.z = w;
		detail_.texture.tex.Create(GL_TEXTURE_3D);
		glTextureStorage3D(detail_.texture.id(), GetMipmapLevels(w, w, w), GL_R8, w, w, w);
		detail_.texture.repeat_size = 5.33f;

		detail_.program = {
			"../shaders/SkyRendering/NoiseGen.comp",
			{{4, 4, 4}, {8, 4, 4}, {8, 8, 4}},
			[](const std::string& src) { return Replace(src, "DETAIL_MAP_GEN", "1"); },
			"DETAIL_MAP"
		};
	}
	{
		const int w = displacement_size;
		displacement_.texture.x = displacement_.texture.y = w;
		displacement_.texture.tex.Create(GL_TEXTURE_2D);
		glTextureStorage2D(displacement_.texture.id(), GetMipmapLevels(w, w), GL_RGBA8, w, w);
		displacement_.texture.repeat_size = 3.51f;

		displacement_.program = {
			"../shaders/SkyRendering/NoiseGen.comp",
			{{16, 8}, {32, 16}, {32, 32}, {8, 8}, {16, 16}, {8, 4}},
			[](const std::string& src) { return Replace(src, "DISPLACEMENT_GEN", "1"); },
			"DISPLACEMENT"
		};
	}

	buffer_.Create();
	glNamedBufferStorage(buffer_.id(), sizeof(VolumetricCloudDefaultMaterialCommonBufferData), NULL, GL_DYNAMIC_STORAGE_BIT);
}

void VolumetricCloudDefaultMaterialCommon::Update(const UpdateParam& param, glm::vec2& additional_delta) {
	cloud_map_.GenerateIfParameterChanged();
	detail_.GenerateIfParameterChanged();
	displacement_.GenerateIfParameterChanged();


	VolumetricCloudDefaultMaterialCommonBufferData buffer;
	buffer.uDensity = density_;

	auto gen_sample_info = [&param, this](const TextureWithInfo& tex, SampleInfo& info, glm::dvec2 offset) {
		info.frequency = 1.0f / tex.repeat_size;
		info.bias = glm::vec2(glm::fract(offset / static_cast<double>(tex.repeat_size)));
		info.k_lod = glm::exp2(lod_bias_) * CalKLod(tex, param.viewport, param.camera);
	};
	const glm::vec2 kLocalWindDirection{ 1.0, 0.0 };
	additional_delta = kLocalWindDirection * wind_speed_ * ImGui::GetIO().DeltaTime;
	auto offset_from_first_cur = param.offset_from_first + glm::dvec2(additional_delta);
	detail_offset_from_first_ += additional_delta * detail_wind_magnify_;
	gen_sample_info(cloud_map_.texture, buffer.uCloudMapSampleInfo, offset_from_first_cur);
	gen_sample_info(detail_.texture, buffer.uDetailSampleInfo, offset_from_first_cur + detail_offset_from_first_);
	gen_sample_info(displacement_.texture, buffer.uDisplacementSampleInfo, offset_from_first_cur);

	glNamedBufferSubData(buffer_.id(), 0, sizeof(buffer), &buffer);
}

void VolumetricCloudDefaultMaterialCommon::Bind() {
	glBindBufferBase(GL_UNIFORM_BUFFER, 3, buffer_.id());

	GLBindTextures({ 
		cloud_map_.texture.id(),
		detail_.texture.id(),
		displacement_.texture.id(),
		}, IVolumetricCloudMaterial::kMaterialTextureUnitBegin);
	GLBindSamplers({ 
		Samplers::Get(Samplers::Wrap::REPEAT, Samplers::Mag::LINEAR, minfilter2d_),
		Samplers::Get(Samplers::Wrap::REPEAT, Samplers::Mag::LINEAR, minfilter3d_),
		Samplers::Get(Samplers::Wrap::REPEAT, Samplers::Mag::LINEAR, minfilter_displacement_),
		}, IVolumetricCloudMaterial::kMaterialTextureUnitBegin);
}

float VolumetricCloudDefaultMaterialCommon::GetSigmaTMax() {
	return density_;
}

void VolumetricCloudDefaultMaterialCommon::DrawGUI() {
	if (ImGui::TreeNode("Cloud Map")) {
		cloud_map_.buffer.DrawGUI();
		ImGui::TreePop();
	}
	if (ImGui::TreeNode("Detail")) {
		detail_.buffer.DrawGUI();
		ImGui::TreePop();
	}
	if (ImGui::TreeNode("Displacement")) {
		displacement_.buffer.DrawGUI();
		ImGui::TreePop();
	}
	ImGui::EnumSelect("Sampler 2D Filter", &minfilter2d_);
	ImGui::EnumSelect("Sampler 3D Filter", &minfilter3d_);
	ImGui::EnumSelect("Sampler Displacement Filter", &minfilter_displacement_);
	ImGui::SliderFloat("Lod Bias", &lod_bias_, -2.0f, 10.0f);
	ImGui::SliderFloat("Cloud Map Repeat Size", &cloud_map_.texture.repeat_size, 0.1f, 100.0f);
	ImGui::SliderFloat("Detail Repeat Size", &detail_.texture.repeat_size, 0.1f, 100.0f);
	ImGui::SliderFloat("Displacement Repeat Size", &displacement_.texture.repeat_size, 0.1f, 100.0f);
	static float density_slider_max = 100.0f;
	ImGui::SliderFloat("Density Slider Max", &density_slider_max, 20.0f, 1000.0f);
	ImGui::SliderFloat("Density", &density_, 0.0f, density_slider_max);
	ImGui::SliderFloat("Wind Speed", &wind_speed_, 0.0f, 1.0f);
	ImGui::SliderFloat("Detail Wind Magnify", &detail_wind_magnify_, -1.0f, 5.0f);
}

void NoiseCreateInfo::DrawGUI() {
	ImGui::InputInt("Seed", &seed);
	ImGui::SliderInt("Base Frequency", &base_frequency, 1, 16);
	ImGui::SliderFloat("Remap Min", &remap_min, 0.0f, 1.0f);
	ImGui::SliderFloat("Remap Max", &remap_max, 0.0f, 1.0f);
}

void CloudMapBuffer::DrawGUI() {
	if (ImGui::TreeNode("Density")) {
		uDensity.DrawGUI();
		ImGui::TreePop();
	}
	if (ImGui::TreeNode("Height")) {
		uHeight.DrawGUI();
		ImGui::TreePop();
	}
}

void DetailBuffer::DrawGUI() {
	if (ImGui::TreeNode("Perlin")) {
		uPerlin.DrawGUI();
		ImGui::TreePop();
	}
	if (ImGui::TreeNode("Worley")) {
		uWorley.DrawGUI();
		ImGui::TreePop();
	}
}

void DisplacementBuffer::DrawGUI() {
	uPerlin.DrawGUI();
}

struct VolumetricCloudDefaultMaterial0BufferData {
	glm::vec2 uDetailParam;
	float uDisplacementScale;
	float padding1;
};

VolumetricCloudDefaultMaterial0::VolumetricCloudDefaultMaterial0() {
	buffer_.Create();
	glNamedBufferStorage(buffer_.id(), sizeof(VolumetricCloudDefaultMaterial0BufferData), NULL, GL_DYNAMIC_STORAGE_BIT);
}

std::string VolumetricCloudDefaultMaterial0::ShaderSrc() {
	return ReadWithPreprocessor("../shaders/SkyRendering/VolumetricCloudDefaultMaterial0.glsl");
}

void VolumetricCloudDefaultMaterial0::Update(const UpdateParam& param, glm::vec2& additional_delta) {
	materail_common_.Update(param, additional_delta);

	VolumetricCloudDefaultMaterial0BufferData buffer;
	buffer.uDetailParam = detail_param_;
	buffer.uDisplacementScale = displacement_scale_;
	glNamedBufferSubData(buffer_.id(), 0, sizeof(buffer), &buffer);
}

void VolumetricCloudDefaultMaterial0::Bind() {
	materail_common_.Bind();
	glBindBufferBase(GL_UNIFORM_BUFFER, 4, buffer_.id());
}

float VolumetricCloudDefaultMaterial0::GetSigmaTMax() {
	return materail_common_.GetSigmaTMax();
}

void VolumetricCloudDefaultMaterial0::DrawGUI() {
	materail_common_.DrawGUI();

	ImGui::SliderFloat("Detail Param X", &detail_param_.x, -1.0f, 1.0f);
	ImGui::SliderFloat("Detail Param Y", &detail_param_.y, 0.0f, 1.0f);
	ImGui::SliderFloat("Displacement Scale", &displacement_scale_, 0.0f, 1.0f);
}

std::string VolumetricCloudDefaultMaterial1::ShaderSrc() {
	return ReadWithPreprocessor("../shaders/SkyRendering/VolumetricCloudDefaultMaterial1.glsl");
}

struct VolumetricCloudDefaultMaterial1BufferData{
	float uBaseDensityThreshold;
	float uBaseHeightHardness;
	float uBaseEdgeHardness;
	float uDetailBase;
	float uDetailScale;
	float uHeightCut;
	float uEdgeCur;
	float padding1;
};

VolumetricCloudDefaultMaterial1::VolumetricCloudDefaultMaterial1() {
	buffer_.Create();
	glNamedBufferStorage(buffer_.id(), sizeof(VolumetricCloudDefaultMaterial1BufferData), NULL, GL_DYNAMIC_STORAGE_BIT);
}

void VolumetricCloudDefaultMaterial1::Update(const UpdateParam& param, glm::vec2& additional_delta) {
	materail_common_.Update(param, additional_delta);

	VolumetricCloudDefaultMaterial1BufferData buffer;
	buffer.uBaseDensityThreshold = base_density_threshold_;
	buffer.uBaseHeightHardness = base_height_hardness_;
	buffer.uBaseEdgeHardness = base_edge_hardness_;
	buffer.uDetailBase = detail_base_;
	buffer.uDetailScale = detail_scale_;
	buffer.uHeightCut = 1.0f - height_cut_;
	buffer.uEdgeCur = edge_cut_;
	glNamedBufferSubData(buffer_.id(), 0, sizeof(buffer), &buffer);
}

void VolumetricCloudDefaultMaterial1::Bind() {
	materail_common_.Bind();
	glBindBufferBase(GL_UNIFORM_BUFFER, 4, buffer_.id());
}

float VolumetricCloudDefaultMaterial1::GetSigmaTMax() {
	return materail_common_.GetSigmaTMax();
}

void VolumetricCloudDefaultMaterial1::DrawGUI() {
	materail_common_.DrawGUI();

	ImGui::SliderFloat("Base Density Threshhold", &base_density_threshold_, 0.0f, 1.0f);
	ImGui::SliderFloat("Base Height Hardness", &base_height_hardness_, 0.0f, 10.0f);
	ImGui::SliderFloat("Base Edge Hardness", &base_edge_hardness_, 0.0f, 10.0f);
	ImGui::SliderFloat("Detail Base", &detail_base_, 0.0f, 1.0f);
	ImGui::SliderFloat("Detail Scale", &detail_scale_, 0.0f, 8.0f);
	ImGui::SliderFloat("Height Cut", &height_cut_, 0.0f, 1.0f);
	ImGui::SliderFloat("Edge Cut", &edge_cut_, 0.0f, 1.0f);
}
