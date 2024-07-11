#include "GLProgram.h"

#include <iostream>
#include <string>
#include <sstream>
#include <array>
#include <stdexcept>
#include <vector>
#include <filesystem>
#include <regex>

#include "Utils.h"

namespace {

std::vector<std::string> SplitLines(const std::string& src, char c) {
	std::vector<std::string> res;
	size_t i = 0;
	size_t j = 0;
	for (;;) {
		j = src.find(c, i);
		if (j == std::string::npos)
			break;
		res.push_back(src.substr(i, j - i));
		i = j + 1;
	}
	res.push_back(src.substr(i, src.size() - 1));
	return res;
}

std::string AnalyzeErrorInfo(const std::string& src, const std::string& info) {
	//return info;
	std::regex pattern(R"(0\((\d+)\) : error C.+?\n)");

	// https://en.cppreference.com/w/cpp/regex/regex_iterator
	auto words_begin = std::sregex_iterator(info.begin(), info.end(), pattern);
	auto words_end = std::sregex_iterator();
	std::stringstream ss;

	auto lines = SplitLines(src, '\n');
	for (auto i = words_begin; i != words_end; ++i) {
		auto match = *i;
		auto match_str = match.str();
		int line = std::stoi(match.str(1));
		auto index = static_cast<size_t>(line) - 1;
		if (index >= lines.size())
			throw std::runtime_error("Unexpected error line number");
		ss << "    " << match_str;
		ss << "    |\n";
		if (index > 0)
			ss << "    | " << lines[index - 1] << "\n";
		ss << "  * | " << lines[index] << "\n";
		if (index + 1 < lines.size())
			ss << "    | " << lines[index + 1] << "\n";
		ss << "    |\n";
	}
	return ss.str();
}

class Shader {
public:
	Shader(GLenum type, const char* src, const std::string& identifier) {
		id_ = glCreateShader(type);
		glShaderSource(id_, 1, &src, nullptr);
		glCompileShader(id_);

		GLint success;
		std::array<char, 4096> info{};
		glGetShaderiv(id_, GL_COMPILE_STATUS, &success);
		if (!success) {
			glGetShaderInfoLog(id_, static_cast<GLsizei>(info.size()), nullptr, info.data());
			glDeleteShader(id_);
			throw std::runtime_error(std::string("Error while compiling ") 
				+ identifier + ":\n" + AnalyzeErrorInfo(src, info.data()));
		}
	}

	~Shader() {
		if (id_ != 0)
			glDeleteShader(id_);
	}

	GLuint Id()const { return id_; }

private:
	GLuint id_ = 0;
};

GLuint CreateProgram(const std::vector<GLuint>& shaders) {
	auto program = glCreateProgram();
	for (auto shader : shaders)
		glAttachShader(program, shader);
	glLinkProgram(program);

	GLint success;
	std::array<char, 1024> info{};
	glGetProgramiv(program, GL_LINK_STATUS, &success);
	if (!success) {
		glGetProgramInfoLog(program, static_cast<GLsizei>(info.size()), NULL, info.data());
		glDeleteProgram(program);
		throw std::runtime_error(std::string("Error while linking program:\n") + info.data());
	}
	return program;
}

}

GLProgram::GLProgram(const char* vertex_src, const char* fragment_src, const Params& params) {
#define STR(NAME) #NAME
	Shader vertex(GL_VERTEX_SHADER, vertex_src, std::string("Vertex Shader: ") + params.tag);
	Shader fragment(GL_FRAGMENT_SHADER, fragment_src, std::string("Fragment Shader: ") + params.tag);
	
	std::vector shaders{ vertex.Id(), fragment.Id() };
	id_ = CreateProgram(shaders);
}

GLProgram::GLProgram(const char* compute_src, const Params& params) {
	Shader compute(GL_COMPUTE_SHADER, compute_src, std::string("Compute Shader: ") + params.tag);

	std::vector shaders{ compute.Id() };
	id_ = CreateProgram(shaders);
}

GLProgram::~GLProgram() {
	if (id_ != 0)
		glDeleteProgram(id_);
}

std::string Replace(std::string src, const std::string& from, const std::string& to) {
	for (;;) {
		auto index = src.find(from);
		if (index == std::string::npos)
			break;
		src.replace(index, from.size(), to);
	}
	return src;
}

std::string ReadWithPreprocessor(const char* filepath) {
	namespace fs = std::filesystem;
	fs::path fpath = filepath;
	auto dir = fpath.parent_path();
	auto src = ReadFile(filepath);
	std::regex pattern("#include[ \t]*\"(.+)\"");
	std::smatch match;
	while (std::regex_search(src, match, pattern)) {
		auto header_path = dir / match.str(1);
		auto header = ReadWithPreprocessor(header_path.string().c_str());
		src = src.substr(0, match.position()) + header 
			+ src.substr(match.position() + match.length(), src.size() - match.position() - match.length());
	}
	return src;
}
