#version 460
TAG_CONF
#include "Common.glsl"
#include "Blur.glsl"
in vec2 vTexCoord;
layout(location = 0) out vec4 FragColor;
layout(binding = 0) uniform sampler2D tex;
layout(binding = 2) uniform sampler2D bloom_mid_tex;
#if DITHER_ENABLE
layout(binding = 3) uniform sampler2D blue_noise;
#endif
layout(location = 0) uniform float filter_width;
layout(location = 1) uniform float bloom_intensity;
layout(location = 2) uniform float exposure;

vec3 ToneMapping(vec3 luminance, float exposure) {
#if TONE_MAPPING == 0
	return 1 - exp(-exposure * luminance);
#elif TONE_MAPPING == 1
	const float k = 10.0 / 16.0; // Make CEToneMapping ACESToneMapping produce similar result with same exposure.

	// https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
    const float A = 2.51f * k * k;
    const float B = 0.03f * k;
    const float C = 2.43f * k * k;
    const float D = 0.59f * k;
    const float E = 0.14f;

    luminance *= exposure;
    return (luminance * (A * luminance + B)) / (luminance * (C * luminance + D) + E);
#endif
}

void main() {
	vec3 luminance = texelFetch(tex, ivec2(gl_FragCoord.xy), 0).rgb;
	vec2 wh = textureSize(tex, 0);
	luminance += Blur(bloom_mid_tex, vTexCoord.xy, vec2(filter_width * (wh.y / wh.x), 0)).rgb * bloom_intensity;
	vec3 color = ToneMapping(luminance, exposure);
#if DITHER_ENABLE
	vec3 srgb_color = Linear2sRGB(color);
	vec3 color_up_bound = sRGB2Linear(srgb_color + 1.0 / 255.0);
	color = mix(color, color_up_bound, texelFetch(blue_noise, ivec2(gl_FragCoord.xy) & 0x3f, 0).x);
#endif
	FragColor = vec4(color, 1.0);
}