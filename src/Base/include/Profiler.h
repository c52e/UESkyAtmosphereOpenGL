#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <optional>

#include "Singleton.h"
#include "gl.hpp"

class Profiler : public Singleton<Profiler> {
public:
	friend Singleton<Profiler>;

	bool profiling = false;

	void NewFrame();

	void Push(std::string label);

	void Pop();

	void Clear();

	void DrawImGui();

private:
	bool frame_profiling = false;

	struct ScopeData {
		ScopeData();

		void Begin();

		void End();

		void DrawImGui(bool is_root = true) const;

		void Save(std::ostream& file) const;
		
		GLuint64 begin_time{};
		double frames{};
		double max_delta_ms = 0;
		double min_delta_ms = 999999;
		double sum_delta_ms = 0;

		double AvgDeltaMs() const {
			return sum_delta_ms / frames;
		}

		void UpdateStatistics(double delta_ms);

		std::vector<std::unordered_map<std::string, ScopeData>::iterator> children_order;
		std::unordered_map<std::string, ScopeData> children;

		static constexpr int kFramesCnt = 3;
		int cur_index = 0;
		GLQueries<kFramesCnt> queries[2];
		bool valid[kFramesCnt]{};
	};

	ScopeData* FindCurrentScope();

	std::vector<std::string> stack_;
	ScopeData root_;
};

#define CONCAT(a, b) CONCAT_INNER(a, b)
#define CONCAT_INNER(a, b) a ## b

#define PROFILE_SCOPE(label) ProfileScope CONCAT(__profile_scope__, __LINE__)(label);
#define PROFILE_SCOPE_DEBUG_GROUP(label) ProfileScopeDebugGroup CONCAT(__profile_scope_debug_group__, __LINE__)(label);

class ProfileScope {
public:
	ProfileScope(const char* label) {
		Profiler::Instance().Push(label);
	}
	~ProfileScope() {
		Profiler::Instance().Pop();
	}
};

class ProfileScopeDebugGroup {
public:
	ProfileScopeDebugGroup(const char* label)
		: scope(label), group(label) {}

	DebugGroup group;
	ProfileScope scope;
};
