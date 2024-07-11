#pragma once

#include <string>

#include <glad/glad.h>

class GLProgram {
public:
	GLProgram() = default;

	struct Params {
		std::string tag = "";
	};
	GLProgram(const char* vertex_src, const char* fragment_src, const Params& params = {});
	GLProgram(const char* compute_src, const Params& params = {});
	~GLProgram();
	GLProgram(const GLProgram&) = delete;
	GLProgram(GLProgram&& rhs) noexcept
		: GLProgram() {
		swap(rhs);
	}

	GLProgram& operator=(GLProgram rhs) noexcept {
		swap(rhs);
		return *this;
	}

	GLuint id() const noexcept { return id_; }

	void swap(GLProgram& rhs) noexcept {
		using std::swap;
		swap(id_, rhs.id_);
	}

private:
	GLuint id_ = 0;
};

constexpr auto kCommonVertexSrc = R"(
#version 460
layout(location = 0) in vec2 aPos;
out vec2 vTexCoord;
void main() {
	vTexCoord = aPos * 0.5 + 0.5;
	gl_Position = vec4(aPos, 0.0, 1.0);
}
)";

std::string Replace(std::string src, const std::string& from, const std::string& to);

std::string ReadWithPreprocessor(const char* filepath);
