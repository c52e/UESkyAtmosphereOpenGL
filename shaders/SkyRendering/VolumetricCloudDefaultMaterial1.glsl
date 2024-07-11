#include "VolumetricCloudDefaultMaterialCommon.glsl"

layout(std140, binding = 4) uniform VolumetricCloudMaterialBufferData{
    float uBaseDensityThreshold;
    float uBaseHeightHardness;
    float uBaseEdgeHardness;
    float uDetailBase;
    float uDetailScale;
    float uHeightCut;
    float uEdgeCur;
    float padding1;
};

void SampleCloud(CloudSampleParams params, inout CloudSampleResult res) {
    vec3 pos = params.pos;
    float height01 = params.height01;
    vec3 camera_pos = params.camera_pos;

    vec4 uvwlod = GetUVWLod(pos, uCloudMapSampleInfo, camera_pos);
    vec2 cloud_type = textureLod(cloud_map, uvwlod.xy, uvwlod.a).rg;
    
    float density = clamp((cloud_type.r - uBaseDensityThreshold) * uBaseEdgeHardness, 0, 1);
    density *= clamp((1 - height01) * uBaseHeightHardness, 0, 1);
    if (density == 0) {
        res.sigma_t = 0;
        return;
    }

    uvwlod = GetUVWLod(pos, uDetailSampleInfo, camera_pos);
    float detail = textureLod(detail_texture, uvwlod.xyz, uvwlod.a).r;
    detail = (detail + uDetailBase) * uDetailScale;
    detail *= max(clamp(height01 - uHeightCut, 0, 1), clamp(uEdgeCur - cloud_type.r, 0, 1));
    density = clamp(density - detail, 0, 1) * uDensity * height01;
    res.sigma_t = density;

    res.attenuation = 0.8;
    res.phase_g = 0.8;
    res.contribution = max(1 - exp(-res.sigma_t * 0.6), min(1 - height01, 0.7));
}