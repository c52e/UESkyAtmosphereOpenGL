#define SMAA_GLSL_4
#define SMAA_PRESET_HIGH

#ifdef VERTEX

#define SMAA_INCLUDE_PS 0
#include "SMAA.hlsl"

layout(location = 0) in vec2 aPos;
out vec2 texcoord;
out vec4 offset;
void main() {
	texcoord = aPos * 0.5 + 0.5;
	SMAANeighborhoodBlendingVS(texcoord, offset);
	gl_Position = vec4(aPos, 0.0, 1.0);
}

#endif


#ifdef FRAGMENT

#include "SMAA.hlsl"

layout(binding = 0) uniform sampler2D colorTex;
layout(binding = 1) uniform sampler2D blendTex;

in vec2 texcoord;
in vec4 offset;
layout(location = 0) out vec4 res;
void main() {
	res = SMAANeighborhoodBlendingPS(texcoord, offset, colorTex, blendTex);
}

#endif
