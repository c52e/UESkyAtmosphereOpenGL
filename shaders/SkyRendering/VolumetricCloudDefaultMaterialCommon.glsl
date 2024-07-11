layout(binding = MATERIAL_TEXTURE_UNIT_BEGIN + 0) uniform sampler2D cloud_map;
layout(binding = MATERIAL_TEXTURE_UNIT_BEGIN + 1) uniform sampler3D detail_texture;
layout(binding = MATERIAL_TEXTURE_UNIT_BEGIN + 2) uniform sampler2D displacement_texture;

struct SampleInfo {
	vec2 bias;
	float frequency;
	float k_lod;
};

layout(std140, binding = 3) uniform VolumetricCloudMaterialCommonBufferData{
	SampleInfo uCloudMapSampleInfo;
	SampleInfo uDetailSampleInfo;
	SampleInfo uDisplacementSampleInfo;
	vec3 padding0;
	float uDensity;
};

vec4 GetUVWLod(vec3 pos, SampleInfo sample_info, vec3 camera_pos) {
	vec3 uvw = pos * sample_info.frequency + vec3(sample_info.bias, 0.0);
	float lod = log2(sample_info.k_lod * distance(pos, camera_pos));
	return vec4(uvw, lod);
}