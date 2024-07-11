layout(binding = MATERIAL_TEXTURE_UNIT_BEGIN + 0) uniform sampler3D voxel;

layout(std140, binding = 3) uniform VolumetricCloudMaterialBufferData{
	vec2 uSampleFrequency;
	float uSampleLodK;
	float uDensity;
	vec2 uSampleBias;
	vec2 voxel_material_padding;
};

void SampleCloud(CloudSampleParams params, inout CloudSampleResult res) {
	vec3 pos = params.pos;
	float height01 = params.height01;
	vec3 camera_pos = params.camera_pos;

    vec2 uv = pos.xy * uSampleFrequency + uSampleBias;
    float lod = log2(uSampleLodK * distance(pos, camera_pos));
    float density = textureLod(voxel, vec3(uv, height01), lod).r;
	res.sigma_t = density * uDensity;
}