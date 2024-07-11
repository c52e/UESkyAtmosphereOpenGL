#include "Profiler.h"

#include <filesystem>
#include <fstream>

#include "imgui.h"

static GLuint64 GetQueryObjectui64(GLuint query) {
	GLint done = 0;
	while (!done)
		glGetQueryObjectiv(query, GL_QUERY_RESULT_AVAILABLE, &done);
	GLuint64 value;
	glGetQueryObjectui64v(query, GL_QUERY_RESULT, &value);
	return value;
}

void Profiler::NewFrame() {
	frame_profiling = profiling;
}

void Profiler::Push(std::string label) {
	if (!frame_profiling)
		return;
	ScopeData* pscope = FindCurrentScope();
	auto itr = pscope->children.find(label);
	if (itr == pscope->children.end()) {
		itr = pscope->children.emplace(label, ScopeData{}).first;
		pscope->children_order.push_back(itr);
	}
	pscope = &itr->second;
	stack_.push_back(std::move(label));
	pscope->Begin();
}

void Profiler::Pop() {
	if (!frame_profiling)
		return;
	ScopeData& scope = *FindCurrentScope();
	scope.End();
	stack_.pop_back();
}

void Profiler::Clear() {
	root_ = {};
	profiling = false;
}

void Profiler::DrawImGui() {
	ImGui::Checkbox("Profiling", &profiling);
	ImGui::SameLine();
	if (ImGui::Button("Clear"))
		Clear();
	if (ImGui::Button("Save..")) {
		ImGui::OpenPopup("Save State");
	}
	if (ImGui::BeginPopupModal("Save State", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) {
		static char path[256] = "profile.csv";
		ImGui::InputText("path", path, std::size(path));
		if (ImGui::Button("OK", ImVec2(120, 0))) {
			namespace fs = std::filesystem;
			std::string folder = "profile/";
			if (!fs::exists(folder))
				fs::create_directory(folder);
			std::ofstream file(folder + path);
			file << std::fixed << std::setprecision(3);
			root_.Save(file);
			ImGui::CloseCurrentPopup();
		}
		ImGui::SameLine();
		if (ImGui::Button("Cancel", ImVec2(120, 0))) {
			ImGui::CloseCurrentPopup();
		}
		ImGui::EndPopup();
	}
	root_.DrawImGui();
}

Profiler::ScopeData* Profiler::FindCurrentScope() {
	ScopeData* res = &root_;
	for (const auto& label : stack_) {
		res = &res->children.at(label);
	}
	return res;
}

Profiler::ScopeData::ScopeData() {
	queries[0].Create(GL_TIMESTAMP);
	queries[1].Create(GL_TIMESTAMP);
}

void Profiler::ScopeData::Begin() {
	if (valid[cur_index])
		begin_time = GetQueryObjectui64(queries[0][cur_index]);
	glQueryCounter(queries[0][cur_index], GL_TIMESTAMP);
}

void Profiler::ScopeData::End() {
	if (valid[cur_index]) {
		auto end = GetQueryObjectui64(queries[1][cur_index]);
		auto elapsed = end - begin_time;
		auto delta_ms = static_cast<double>(elapsed) * 1e-6;
		UpdateStatistics(delta_ms);
	}
	glQueryCounter(queries[1][cur_index], GL_TIMESTAMP);
	valid[cur_index] = true;
	cur_index = (cur_index + 1) % kFramesCnt;
}

void Profiler::ScopeData::DrawImGui(bool is_root) const {
	if (is_root) {
		ImGui::Text("avg\t\tmax\t\tmin"
			, AvgDeltaMs(), max_delta_ms, min_delta_ms);
	}
	else {
		ImGui::Text("%.3f ms\t%.3f ms\t%.3f ms"
			, AvgDeltaMs(), max_delta_ms, min_delta_ms);
	}
	double children_sum = 0;
	for (const auto& itr : children_order) {
		const auto& [label, scope] = *itr;
		children_sum += scope.AvgDeltaMs();
		if (scope.children.empty()) {
			ImGui::Text("%.3f ms\t%.3f ms\t%.3f ms (%s)"
				, scope.AvgDeltaMs(), scope.max_delta_ms, scope.min_delta_ms, label.c_str());
		}
		else if (ImGui::TreeNode(label.c_str())) {
			scope.DrawImGui(false);
			ImGui::TreePop();
		}
	}
	ImGui::Text("%.3f ms\t(children sum)", children_sum);
}

void Profiler::ScopeData::Save(std::ostream& file) const {
	for (const auto& itr : children_order) {
		const auto& [label, scope] = *itr;
		if (scope.children.empty())
			file << label << "," << scope.AvgDeltaMs() << "\n";
		else
			scope.Save(file);
	}
}

void Profiler::ScopeData::UpdateStatistics(double delta_ms) {
	frames++;
	min_delta_ms = std::min(min_delta_ms, delta_ms);
	max_delta_ms = std::max(max_delta_ms, delta_ms);
	sum_delta_ms += delta_ms;
}
