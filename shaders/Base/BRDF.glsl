#include "Common.glsl"

#ifdef SPECTRUM_RENDERING
#define SpectrumType float
#else
#define SpectrumType vec3
#endif

struct MaterialData {
	SpectrumType F0;
	SpectrumType diffuse;
	vec3 normal;
	float roughness;
	float ao;
};

MaterialData LoadMeterialData(sampler2D albedo_tex, sampler2D normal_tex, sampler2D orm_tex, ivec2 xy) {
	SpectrumType albedo = SpectrumType(texelFetch(albedo_tex, xy, 0));
	vec3 normal = texelFetch(normal_tex, xy, 0).rgb;
	vec3 orm = texelFetch(orm_tex, xy, 0).rgb;

	MaterialData data;
	float metallic = orm.b;
	data.F0 = mix(SpectrumType(0.04), albedo, metallic);
	data.diffuse = albedo - albedo * metallic;
	data.normal = normal;
	data.roughness = orm.g;
	data.ao = orm.r;

	return data;
}

//--------------------------------------
// from UE4
float Pow5(float x) {
	float x2 = x * x;
	return x2 * x2 * x;
}

SpectrumType F_Schlick(float HdotV, SpectrumType F0) {
	return F0 + (1.0 - F0) * Pow5(1.0 - HdotV);
}

float D_GGX(float a, float NdotH) {
	float a2 = a * a;
	float d = (NdotH * a2 - NdotH) * NdotH + 1.0;
	return a2 / (PI * d * d);
}

float Vis_SmithJointApprox(float a, float NdotV, float NdotL) {
	float Vis_SmithV = NdotL * (NdotV * (1 - a) + a);
	float Vis_SmithL = NdotV * (NdotL * (1 - a) + a);
	return 0.5 / max(Vis_SmithV + Vis_SmithL, 1e-9);
}

vec3 ImportanceSampleGGX(vec2 E, float a) {
	float a2 = a * a;
	float Phi = 2 * PI * E.x;
	float CosTheta = sqrt((1 - E.y) / (1 + (a2 - 1) * E.y));
	float SinTheta = sqrt(1 - CosTheta * CosTheta);

	vec3 H;
	H.x = SinTheta * cos(Phi);
	H.y = SinTheta * sin(Phi);
	H.z = CosTheta;

	return H;
}

vec3 ImportanceSampleGGX(vec2 E, float a, vec3 N) {
	vec3 H = ImportanceSampleGGX(E, a);
	vec3 t0, t1;
	CreateOrthonormalBasis(N, t0, t1);
	return t0 * H.x + t1 * H.y + N * H.z;
}
//--------------------------------------

SpectrumType BRDF(float NdotL, float NdotV, float NdotH, float HdotV, MaterialData data) {
	float a = data.roughness * data.roughness;
	float D = min(D_GGX(a, NdotH), 1e9);
	float Vis = Vis_SmithJointApprox(a, NdotV, NdotL);
	SpectrumType F = F_Schlick(HdotV, data.F0);

	SpectrumType specular = SpectrumType(D * Vis);
	SpectrumType diffuse = INV_PI * data.diffuse;
	return mix(diffuse, specular, F);
}

// https://cseweb.ucsd.edu/~ravir/papers/envmap/envmap.pdf
vec3 GetSHIrradiance(vec3 N, vec4 Llm[9]) {
	const float c1 = 0.429043;
	const float c2 = 0.511664;
	const float c3 = 0.743125;
	const float c4 = 0.886227;
	const float c5 = 0.247708;

	vec3 L00 = Llm[0].rgb;
	vec3 L1_1 = Llm[1].rgb;
	vec3 L10 = Llm[2].rgb;
	vec3 L11 = Llm[3].rgb;
	vec3 L2_2 = Llm[4].rgb;
	vec3 L2_1 = Llm[5].rgb;
	vec3 L20 = Llm[6].rgb;
	vec3 L21 = Llm[7].rgb;
	vec3 L22 = Llm[8].rgb;

	float x = N.x;
	float y = N.y;
	float z = N.z;

	return c1 * (x * x - y * y) * L22 + c3 * (z * z) * L20 + c4 * L00 - c5 * L20
		+ 2.0 * c1 * (x * y * L2_2 + x * z * L21 + y * z * L2_1)
		+ 2.0 * c2 * (x * L11 + y * L1_1 + z * L10);
}

#ifndef SPECTRUM_RENDERING
vec3 GetAmbient(sampler2D env_brdf_lut, float roughness_lod_max, samplerCube prefiltered_radiance_texture,
		vec3 V, MaterialData data, vec4 Llm[9]) {
	vec3 N = data.normal;
	float NdotV = clamp(dot(N, V), 0, 1);
#if 0
	vec3 F = F_Schlick(NdotV, data.F0);
#else
	vec3 F = F_Schlick(NdotV * 0.8 + 0.2, data.F0);
#endif
#if 0
	vec3 approx_irradiance_over_pi = textureLod(prefiltered_radiance_texture, N, roughness_lod_max).rgb;
#else
	vec3 approx_irradiance_over_pi = GetSHIrradiance(N, Llm) * INV_PI;
#endif
	vec3 diffuse = (data.diffuse - data.diffuse * F) * approx_irradiance_over_pi;

	vec3 R = 2.0 * NdotV * N - V;
	vec3 prefiltered_radiance = textureLod(prefiltered_radiance_texture, R, data.roughness * roughness_lod_max).rgb;
	vec2 brdf = texture(env_brdf_lut, vec2(NdotV, data.roughness)).rg;
	vec3 specular = prefiltered_radiance * (data.F0 * brdf.x + brdf.y);

	return data.ao * (diffuse + specular);
}
#endif

vec3 CosineSample(vec2 E, vec3 N, out float pdf) {
#if 0
	float sin_theta = sqrt(E.x);
	float cos_theta = sqrt(clamp(1.0 - sin_theta * sin_theta, 0, 1));
	vec3 t0, t1;
	CreateOrthonormalBasis(N, t0, t1);
	float phi = 2.0 * PI * E.y;
	pdf = INV_PI * cos_theta;
	return vec3(sin_theta * sin(phi) * t0 + sin_theta * cos(phi) * t1 + cos_theta * N);
#else
	// https://cseweb.ucsd.edu/~tzli/cse272/wi2023/lectures/28_misc.pdf
	float theta = 2.0 * PI * E.x;
	E.y = 2.0 * E.y - 1.0;
	vec3 spherePoint = vec3(sqrt(1.0 - E.y * E.y) * vec2(cos(theta), sin(theta)), E.y);
	vec3 res = normalize(N + spherePoint);
	pdf = dot(res, N);
	return res;
#endif
}
