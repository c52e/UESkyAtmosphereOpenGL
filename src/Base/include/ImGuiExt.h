#pragma once

#include <magic_enum.hpp>
#include <imgui.h>

namespace ImGui {

    namespace {
        template<class T>
        constexpr size_t MaxEnumStringViewSize() {
            size_t maxsize = 0;
            for (const auto& strview : magic_enum::enum_names<T>()) {
                maxsize = maxsize < strview.size() ? strview.size() : maxsize;
            }
            return maxsize;
        }
    }

    template<class Enum, class=std::enable_if_t<std::is_enum_v<Enum> && (magic_enum::enum_count<Enum>() > 0)>>
    void EnumSelect(const char* name, Enum* penum) {
        static std::array<std::string, magic_enum::enum_count<Enum>()> enums;
        if (enums[0].empty()) {
            constexpr auto& names = magic_enum::enum_names<Enum>();
            auto i = enums.begin();
            auto j = names.begin();
            while (i < enums.end())
                //*i++ = *j++; // doesn't work in msvc c++20 lastest
                *i++ = std::string(*j++);
        }
        char buf[MaxEnumStringViewSize<Enum>() + 1];
        auto curstrview = magic_enum::enum_name(*penum);
        memcpy(buf, &curstrview[0], curstrview.size());
        buf[curstrview.size()] = 0;
        if (ImGui::BeginCombo(name, buf)) {
            for (const auto& enumstr : enums) {
                bool is_selected = strcmp(buf, enumstr.c_str()) == 0;
                if (ImGui::Selectable(enumstr.c_str(), is_selected))
                    *penum = magic_enum::enum_cast<Enum>(enumstr).value();
                if (is_selected)
                    ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
        }
    }

}
