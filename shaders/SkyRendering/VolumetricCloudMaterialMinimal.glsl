layout(std140, binding = 3) uniform VolumetricCloudMaterialBufferData {
	vec3 padding;
	float uDensity;
};

void SampleCloud(CloudSampleParams params, inout CloudSampleResult res) {
    res.sigma_t = uDensity;
}