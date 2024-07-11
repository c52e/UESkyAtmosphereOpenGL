#define _CRT_SECURE_NO_WARNINGS
#define TINYGLTF_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <tiny_gltf.h>

namespace {
	struct StbiImageFlip{
		StbiImageFlip() {
			stbi_set_flip_vertically_on_load(true);
			stbi_flip_vertically_on_write(true);
		}
	} _;
}
