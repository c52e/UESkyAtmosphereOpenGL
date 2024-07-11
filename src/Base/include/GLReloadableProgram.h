#pragma once

#include <functional>
#include <memory>
#include <type_traits>

#include <glm/glm.hpp>

#include "GLProgram.h"
#include "ObjectsSet.h"

class GLReloadableProgram: private ObjectsSet<GLReloadableProgram> {
public:
	friend ObjectsSet<GLReloadableProgram>;

	GLReloadableProgram() = default;

	GLReloadableProgram(GLReloadableProgram&&) = default;

	template<class T, class=std::enable_if_t<std::is_constructible_v<std::function<GLProgram()>, T>>>
	GLReloadableProgram(T&& func) : loader_(std::forward<T>(func)) {
		Reload();
	}

	GLReloadableProgram& operator=(GLReloadableProgram&&) = default;

	GLuint id() const noexcept { return program_.id(); }

	void Reload();

	static void ReloadAll();
	
private:
	std::function<GLProgram()> loader_;
	GLProgram program_;
};

class GLReloadableComputeProgram : private ObjectsSet<GLReloadableComputeProgram> {
public:
	friend ObjectsSet<GLReloadableComputeProgram>;

	GLReloadableComputeProgram() = default;

	GLReloadableComputeProgram(std::string path, const std::vector<glm::ivec2>& localsizes
			, std::function<std::string(const std::string&)> post_process = [](const std::string& src) { return src; }
			, std::string tag = "")
		: GLReloadableComputeProgram(std::move(path), ToVec3(localsizes), std::move(post_process), std::move(tag)) {}

	GLReloadableComputeProgram(std::string path, std::vector<glm::ivec3> localsizes
			, std::function<std::string(const std::string&)> post_process = [](const std::string& src) { return src; }
			, std::string tag = "");

	GLReloadableComputeProgram(GLReloadableComputeProgram&&) = default;

	GLReloadableComputeProgram& operator=(GLReloadableComputeProgram&&) = default;

	void Dispatch(const glm::ivec2& globalsize) const {
		glm::ivec3 globalsize3(globalsize, 1);
		Dispatch(globalsize3);
	}

	void Dispatch(const glm::ivec3& globalsize) const {
		const auto& localsize = data_->localsizes[data_->index];
		auto groupsize = (globalsize + localsize - 1) / localsize;
		glDispatchCompute(groupsize.x, groupsize.y, groupsize.z);
	}

	GLuint id() const { return data_->programs[data_->index].id(); }

	void DrawGUI();

	static void DrawGUIAll();

private:
	void Construct(int i);

	static std::vector<glm::ivec3> ToVec3(const std::vector<glm::ivec2>& in);

	struct Data {
		std::string path;
		std::string display_text;
		std::function<std::string(const std::string&)> post_process;
		std::vector<GLReloadableProgram> programs;
		std::vector<glm::ivec3> localsizes;
		std::vector<std::string> localsizes_str;
		int index = 0;
	};
	std::unique_ptr<Data> data_;
};
