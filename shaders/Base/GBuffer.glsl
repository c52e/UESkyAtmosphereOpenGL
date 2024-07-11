#version 460

#if VERTEX
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aNormal;
layout(location = 2) in vec2 aUv;
layout(location = 3) in vec3 aTangent;
out mat3 vNormalMatrix;
out vec2 vUv;
layout(location = 0) uniform mat4 mvp;
layout(location = 1) uniform mat3 normal_matrix;
void main() {
	gl_Position = mvp * vec4(aPos, 1.0);
	vec3 N = aNormal;
	vec3 T = aTangent;

	// Use glTF (OpenGL) normal direction
	// https://github.com/KhronosGroup/glTF/issues/952
	vec3 B = cross(N, T);
	vNormalMatrix = normal_matrix * mat3(T, B, N);
	vUv = aUv;
}
#endif

#if FRAGMENT
in mat3 vNormalMatrix;
in vec2 vUv;
layout(location = 0) out vec4 Albedo;
layout(location = 1) out vec4 Normal;
layout(location = 2) out vec4 ORM;
layout(location = 2) uniform vec3 albedo_factor;
layout(location = 3) uniform vec2 metallic_roughness_factor;
layout(location = 4) uniform float ao_factor;
layout(binding = 0) uniform sampler2D albedo_texture;
layout(binding = 1) uniform sampler2D normal_texture;
layout(binding = 2) uniform sampler2D metallic_roughness_texture;
layout(binding = 3) uniform sampler2D ao_texture;

void main() {
	vec3 albedo = albedo_factor * texture(albedo_texture, vUv).rgb;
	vec2 metallic_roughness = metallic_roughness_factor * texture(metallic_roughness_texture, vUv).rg;
	float ao = ao_factor * texture(ao_texture, vUv).r;
	Albedo = vec4(albedo, 1.0);
	vec3 tangent_normal = texture(normal_texture, vUv).xyz * 2.0 - 1.0;
	vec3 normal = normalize(vNormalMatrix * tangent_normal);
	Normal = vec4(normal, 1.0);
	ORM = vec4(ao, metallic_roughness.gr, 1.0);
}
#endif
