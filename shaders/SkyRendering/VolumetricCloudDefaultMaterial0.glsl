#include "VolumetricCloudDefaultMaterialCommon.glsl"

layout(std140, binding = 4) uniform VolumetricCloudMaterialBufferData{
    vec2 uDetailParam;
    float uDisplacementScale;
    float padding1;
};

float CalHeightMask(float cloud_type, float height01) {
    float height_in_type = clamp(height01 / cloud_type, 0, 1);
    return clamp(height_in_type * (height_in_type - 1.0) * -4.0, 0, 1);
}

float Remap01(float x, float x0, float x1) {
    return clamp((x - x0) / (x1 - x0), 0, 1);
}

void SampleCloud(CloudSampleParams params, inout CloudSampleResult res) {
    vec3 pos = params.pos;
    float height01 = params.height01;
    vec3 camera_pos = params.camera_pos;

    vec4 uvwlod = GetUVWLod(pos, uCloudMapSampleInfo, camera_pos);
    vec2 cloud_type = textureLod(cloud_map, uvwlod.xy, uvwlod.a).rg;
    
    vec3 displace_vector = vec3(0);
    uvwlod = GetUVWLod(pos, uDisplacementSampleInfo, camera_pos);
    displace_vector.xy += textureLod(displacement_texture, uvwlod.xy, uvwlod.a).rg;
    displace_vector.xz += textureLod(displacement_texture, uvwlod.xz, uvwlod.a).ba;

    pos += uDisplacementScale * displace_vector;
    uvwlod = GetUVWLod(pos, uDetailSampleInfo, camera_pos);
    float detail = textureLod(detail_texture, uvwlod.xyz, uvwlod.a).r;
    detail = detail * uDetailParam.x + uDetailParam.y;
    res.sigma_t = Remap01(cloud_type.r * CalHeightMask(cloud_type.g, height01), detail, 1.0) * height01 * uDensity;

    res.attenuation = 0.8;
    res.phase_g = 0.8;
    res.contribution = max(1 - exp(-res.sigma_t * 0.8), min(1 - height01, 0.7));
}