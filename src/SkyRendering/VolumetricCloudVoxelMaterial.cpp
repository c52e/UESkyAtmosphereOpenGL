#include "VolumetricCloudVoxelMaterial.h"

#if __has_include(<openvdb/openvdb.h>)
#define HAS_INCLUDE_OPENVDB 1
#include <openvdb/openvdb.h>
#else
#define HAS_INCLUDE_OPENVDB 0
#include <Windows.h>
#undef max
#undef min
#endif

#include <imgui.h>

#include "ImageLoader.h"
#include "GLProgram.h"
#include "log.h"

struct VolumetricCloudVoxelMaterial::BufferData {
	glm::vec2 uSampleFrequency;
	float uSampleLodK;
	float uDensity;
	glm::vec2 uSampleBias;
	glm::vec2 voxel_material_padding;
};

VolumetricCloudVoxelMaterial::VolumetricCloudVoxelMaterial() {
	buffer_.Create();
	glNamedBufferStorage(buffer_.id(), sizeof(BufferData), NULL, GL_DYNAMIC_STORAGE_BIT);

	sampler_.Create();
	glSamplerParameteri(sampler_.id(), GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glSamplerParameteri(sampler_.id(), GL_TEXTURE_MIN_FILTER, GL_NEAREST_MIPMAP_NEAREST);
	glSamplerParameteri(sampler_.id(), GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
	glSamplerParameteri(sampler_.id(), GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
	glSamplerParameteri(sampler_.id(), GL_TEXTURE_WRAP_R, GL_CLAMP_TO_BORDER);
	float border_color[] = { 0.0f, 0.0f, 0.0f, 0.0f };
	glSamplerParameterfv(sampler_.id(), GL_TEXTURE_BORDER_COLOR, border_color);

#if HAS_INCLUDE_OPENVDB
    openvdb::initialize();
    openvdb::io::File file("../data/wdas/wdas_cloud_quarter.vdb");
	try {
		file.open();
	} catch(openvdb::IoError& e) {
		LOG_INFO("\"wdas_cloud_quarter.vdb\" not found, use \"wdas_cloud_sixteenth.vdb\".")
		file = openvdb::io::File("../data/wdas/wdas_cloud_sixteenth.vdb");
		file.open();
	}
    auto base_grid = file.readGrid(file.beginName().gridName());
    file.close();
    auto grid = openvdb::gridPtrCast<openvdb::FloatGrid>(base_grid);
    auto mincoord = openvdb::Coord::max();
    auto maxcoord = openvdb::Coord::min();
    for (auto iter = grid->beginValueOn(); iter; ++iter) {
        mincoord.minComponent(iter.getBoundingBox().min());
        maxcoord.maxComponent(iter.getBoundingBox().max());
    }
    auto dim = maxcoord - mincoord + openvdb::Coord(1);
    voxel_dim_ = { dim.x(), dim.z(), dim.y() }; // swap yz
    std::vector<float> data(size_t(voxel_dim_.x) * voxel_dim_.y * voxel_dim_.z, 0.0f);

    for (auto iter = grid->beginValueOn(); iter; ++iter) {
		auto value = iter.getValue();
		auto minxyz = iter.getBoundingBox().min();
		auto maxxyz = iter.getBoundingBox().max();
		for (auto x = minxyz.x(); x <= maxxyz.x(); ++x) {
			for (auto y = minxyz.y(); y <= maxxyz.y(); ++y) {
				for (auto z = minxyz.z(); z <= maxxyz.z(); ++z) {
					auto offset = decltype(mincoord){ x,y,z } - mincoord;
					data[size_t(voxel_dim_.x) * voxel_dim_.y * offset.y() + size_t(voxel_dim_.x) * offset.z() + offset.x()] = value;
				}
			}
		}
    }

	voxel_.Create(GL_TEXTURE_3D);
	glTextureStorage3D(voxel_.id(), GetMipmapLevels(voxel_dim_.x, voxel_dim_.y, voxel_dim_.z)
		, GL_R8, voxel_dim_.x, voxel_dim_.y, voxel_dim_.z);
	glTextureSubImage3D(voxel_.id(), 0, 0, 0, 0, voxel_dim_.x, voxel_dim_.y, voxel_dim_.z, GL_RED, GL_FLOAT, data.data());
	glGenerateTextureMipmap(voxel_.id());
#else
#define OPENVDB_NOT_FOUND_MSG \
	"  To make voxel material available, you need to install OpenVDB and apply user-wide integration:\n\n" \
	"    vcpkg install openvdb:x64-windows\n" \
	"    vcpkg integrate install\n\n" \
	"  Then restart Visual Studio and rebuild the project\n\n" \
	"  You can download vcpkg from https://github.com/microsoft/vcpkg \n"
#pragma message(OPENVDB_NOT_FOUND_MSG)
	MessageBoxA(0, OPENVDB_NOT_FOUND_MSG, "Error", MB_OK);
	throw std::runtime_error(OPENVDB_NOT_FOUND_MSG);
#endif
}

std::string VolumetricCloudVoxelMaterial::ShaderSrc() {
	return ReadWithPreprocessor("../shaders/SkyRendering/VolumetricCloudMaterialVoxel.glsl");
}

void VolumetricCloudVoxelMaterial::Update(const UpdateParam& param, glm::vec2& additional_delta) {
	BufferData buffer;
	buffer.uDensity = density_;
	buffer.uSampleFrequency = 1.0f / width_;
	buffer.uSampleBias = glm::vec2((param.offset_from_first + glm::dvec2(base_)) / glm::dvec2(width_));

	auto max_width = static_cast<float>(glm::max(voxel_dim_.x, glm::max(voxel_dim_.y, voxel_dim_.z)));
	buffer.uSampleLodK = glm::exp2(lod_bias_) * max_width * param.camera.tangent_half_fovy() * 2.0f
		/ (std::min(width_.x, width_.y) * static_cast<float>(param.viewport.y));

	glNamedBufferSubData(buffer_.id(), 0, sizeof(buffer), &buffer);
}

void VolumetricCloudVoxelMaterial::Bind() {
	glBindBufferBase(GL_UNIFORM_BUFFER, 3, buffer_.id());

	GLBindTextures({voxel_.id()}, IVolumetricCloudMaterial::kMaterialTextureUnitBegin);
	GLBindSamplers({ sampler_.id() }, IVolumetricCloudMaterial::kMaterialTextureUnitBegin);
}

float VolumetricCloudVoxelMaterial::GetSigmaTMax() {
	return density_;
}

void VolumetricCloudVoxelMaterial::DrawGUI() {
	static float density_slider_max = 100.0f;
	// http://www.specinc.com/extinctiometer-theory
	// Cumulus clouds may have extinction coefficients on the order of 200 - 300 km^-1,
	// while thin cirrus clouds may have an extinction coefficient on the order of 0.5 - 5 km^-1.
	ImGui::SliderFloat("Density Slider Max", &density_slider_max, 20.0f, 1000.0f);
	ImGui::SliderFloat("Density", &density_, 0.0f, density_slider_max);
	ImGui::SliderFloat("Lod Bias", &lod_bias_, -2.0f, 10.0f);
	ImGui::SliderFloat("Base X", &base_.x, -10.0f, 10.0f);
	ImGui::SliderFloat("Base Y", &base_.y, -10.0f, 10.0f);
	ImGui::SliderFloat("Width X", &width_.x, 0.1f, 10.0f);
	ImGui::SliderFloat("Width Y", &width_.y, 0.1f, 10.0f);
}
