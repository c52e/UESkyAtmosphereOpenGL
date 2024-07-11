#pragma once

#include "reflection.h"

#include <memory>
#include <array>
#include <vector>
#include <list>
#include <map>
#include <unordered_map>

#include <imgui.h>
#include <magic_enum.hpp>

#include "util.h"

namespace reflection {

class IAutoImGui : public IReflectionBase<IAutoImGui> {
};

template<>
class IType<IAutoImGui> {
public:
    virtual void DrawAutoImGui(void* addr, const char* name, const UserdataBase* userdata) const = 0;
};

template<class T>
void DrawAutoImGui(T& object, const char* name = nullptr) {
    const Type<IAutoImGui, T>::Userdata userdata = {};
    Type<IAutoImGui, T>::GetIType()->DrawAutoImGui(&object, name, &userdata);
}

template<>
class Type<IAutoImGui, int> : public TypeBase<IAutoImGui, int> {
public:
    struct Userdata : UserdataBase {
        int Min = 0;
        int Max = 10;
    };

    void DrawAutoImGui(void* addr, const char* name, const UserdataBase* userdata) const override {
        auto p = static_cast<ValueType*>(addr);
        auto arg = static_cast<const Userdata*>(userdata);
        ImGui::SliderInt(name, p, arg->Min, arg->Max);
    }
};

template<>
class Type<IAutoImGui, bool> : public TypeBase<IAutoImGui, bool> {
public:
    void DrawAutoImGui(void* addr, const char* name, const UserdataBase* userdata) const override {
        auto p = static_cast<ValueType*>(addr);
        ImGui::Checkbox(name, p);
    }
};

template<>
class Type<IAutoImGui, float> : public TypeBase<IAutoImGui, float> {
public:
    struct Userdata : UserdataBase {
        float Min = 0.0f;
        float Max = 1.0f;
    };

    void DrawAutoImGui(void* addr, const char* name, const UserdataBase* userdata) const override {
        auto p = static_cast<ValueType*>(addr);
        auto arg = static_cast<const Userdata*>(userdata);
        ImGui::SliderFloat(name, p, arg->Min, arg->Max);
    }
};

template<class T>
class Type<IAutoImGui, T, std::enable_if_t<std::is_enum_v<T>>> : public TypeBase<IAutoImGui, T> {
public:
    Type() {
        constexpr auto& names = magic_enum::enum_names<T>();
        auto i = enums.begin();
        auto j = names.begin();
        while (i < enums.end())
            *i++ = *j++;
    }

    void DrawAutoImGui(void* addr, const char* name, const UserdataBase* userdata) const override {
        auto p = static_cast<T*>(addr);
        char buf[MaxEnumStringViewSize<T>() + 1];
        auto curstrview = magic_enum::enum_name(*p);
        memcpy(buf, &curstrview[0], curstrview.size());
        buf[curstrview.size()] = 0;
        if (ScopeImGuiCombo combo(name, buf); combo) {
            for (const auto& enumstr : enums) {
                bool is_selected = strcmp(buf, enumstr.c_str()) == 0;
                if (ImGui::Selectable(enumstr.c_str(), is_selected))
                    *p = magic_enum::enum_cast<T>(enumstr).value();
                if (is_selected)
                    ImGui::SetItemDefaultFocus();
            }
        }
    }
private:
    std::array<std::string, magic_enum::enum_count<T>()> enums;
};

template<class T>
class Type<IAutoImGui, T, std::enable_if_t<std::is_base_of_v<IAutoImGui, T> || IsReflectableStruct<IAutoImGui, T>::value>>
    : public TypeBase<IAutoImGui, T> {
public:
    using ValueType = T;

    void DrawAutoImGui(void* addr, const char* name, const UserdataBase* userdata) const override {
        auto& v = *static_cast<ValueType*>(addr);
        if (ScopeImGuiTreeNode tree(name); tree) {
            for (const auto& [name, fun] : GetFieldTable(v, static_cast<IAutoImGui*>(nullptr))) {
                auto info = fun(&v);
                info.type->DrawAutoImGui(info.address, name.c_str(), info.userdata);
            }
        }
    }
};

template<class _Ty, class _Dx>
class Type<IAutoImGui, std::unique_ptr<_Ty, _Dx>> : public TypeBase<IAutoImGui, std::unique_ptr<_Ty, _Dx>> {
public:
    using ValueType = std::unique_ptr<_Ty, _Dx>;

    struct Userdata : Type< IAutoImGui, _Ty>::Userdata {};

    void DrawAutoImGui(void* addr, const char* name, const UserdataBase* userdata) const override {
        auto& v = *static_cast<ValueType*>(addr);

        if (v) {
            ScopeImGuiTreeNode tree(name);
            if (ScopeImGuiPopupContextItem popup; popup) {
                if (ImGui::MenuItem("delete")) {
                    v.reset();
                }
            }
            if (tree && v != nullptr) {
                Type<IAutoImGui, _Ty>::GetIType()->DrawAutoImGui(v.get(), "value", userdata);
            }
        }
        else {
            ImGui::Text("%s is null", name);
            if (ScopeImGuiPopupContextItem popup(name); popup) {
                if constexpr (std::is_base_of_v<IAutoImGui, _Ty> && SubclassInfo<_Ty>::has) {
                    for (const auto& [name, factory] : SubclassInfo<_Ty>::GetFactoryTable()) {
                        if (ImGui::MenuItem(name.c_str())) {
                            v.reset(factory());
                            break;
                        }
                    }
                }
                else {
                    if (ImGui::MenuItem("new")) {
                        v.reset(new _Ty());
                    }
                }
            }
        }
    }
};

template<class _Ty, size_t _Size>
struct _AutoImGuiArrayTypeHelper {
    static void DrawAutoImGui(void* addr, const char* name, const UserdataBase* userdata) {
        auto arr = static_cast<_Ty*>(addr);

        if (ScopeImGuiTreeNode tree(name);  tree) {
            constexpr size_t kBufSize = 128;
            char buf[kBufSize];
            for (size_t i = 0; i < _Size; ++i) {
                snprintf(buf, kBufSize, "%zu", i);
                Type<IAutoImGui, _Ty>::GetIType()->DrawAutoImGui(&arr[i], buf, userdata);
            }
        }
    }
};

template<class _Ty, size_t _Size>
class Type< IAutoImGui, _Ty[_Size]>
    : public TypeBase<IAutoImGui, _Ty[_Size]> {
public:
    struct Userdata : Type< IAutoImGui, _Ty>::Userdata {};

    void DrawAutoImGui(void* addr, const char* name, const UserdataBase* userdata) const override {
        _AutoImGuiArrayTypeHelper<_Ty, _Size>::DrawAutoImGui(addr, name, userdata);
    }
};

template<class _Ty, size_t _Size>
class Type< IAutoImGui, std::array<_Ty, _Size>>
    : public TypeBase<IAutoImGui, std::array<_Ty, _Size>> {
public:
    struct Userdata : Type< IAutoImGui, _Ty>::Userdata {};

    void DrawAutoImGui(void* addr, const char* name, const UserdataBase* userdata) const override {
        _AutoImGuiArrayTypeHelper<_Ty, _Size>::DrawAutoImGui(addr, name, userdata);
    }
};

template<template<class _Ty, class _Alloc> class ContainerType, class _Ty, class _Alloc>
class Type<IAutoImGui, ContainerType<_Ty, _Alloc>
    , std::enable_if_t<std::is_same_v<ContainerType<_Ty, _Alloc>, std::vector<_Ty, _Alloc>>
                    || std::is_same_v<ContainerType<_Ty, _Alloc>, std::list<_Ty, _Alloc>>>>
    : public TypeBase<IAutoImGui, ContainerType<_Ty, _Alloc>> {
public:
    using ValueType = ContainerType<_Ty, _Alloc>;

    struct Userdata : Type< IAutoImGui, _Ty>::Userdata {};

    void DrawAutoImGui(void* addr, const char* name, const UserdataBase* userdata) const override {
        auto& v = *static_cast<ValueType*>(addr);

        ScopeImGuiTreeNode tree(name);
        if (ScopeImGuiPopupContextItem popup; popup) {
            if (ImGui::MenuItem("clear")) {
                v.clear();
            }
            else if (ImGui::MenuItem("append")) {
                v.emplace_back();
            }
            else if (ImGui::MenuItem("pop")) {
                if (!v.empty()) {
                    v.pop_back();
                }
            }
        }
        if (tree && !v.empty()) {
            constexpr size_t kBufSize = 128;
            char buf[kBufSize];
            size_t i = 0;
            snprintf(buf, kBufSize, "%zu", i++);
            Type<IAutoImGui, _Ty>::GetIType()->DrawAutoImGui(&v.front(), buf, userdata);
            for (auto itr = ++v.begin(); itr != v.end(); ++itr) {
                snprintf(buf, kBufSize, "exchange(%zu, %zu)", i - 1, i);
                if (ImGui::Button(buf)) {
                    auto pre = itr; --pre;
                    std::swap(*itr, *pre);
                    break;
                }
                snprintf(buf, kBufSize, "%zu", i++);
                Type<IAutoImGui, _Ty>::GetIType()->DrawAutoImGui(&*itr, buf, userdata);
            }
        }
    }
};

template<class T>
struct _AutoImGuiMapTypeHelper {
    using _Ty = typename T::mapped_type;

    static void DrawAutoImGui(void* addr, const char* name, const UserdataBase* userdata) {
        auto& v = *static_cast<T*>(addr);

        ScopeImGuiTreeNode tree(name);
        if (ScopeImGuiPopupContextItem popup; popup) {
            static char keybuf[128];
            if (ImGui::MenuItem("clear")) {
                v.clear();
            }
            ImGui::InputText("##keyinput", keybuf, IM_ARRAYSIZE(keybuf));
            ImGui::SameLine();
            if (ImGui::Button("add")) {
                std::pair<const char*, _Ty> pair;
                pair.first = keybuf;
                v.insert(std::move(pair));
            }
            
        }
        if (tree) {
            int id = 0;
            for (auto itr = v.begin(); itr != v.end(); ++itr) {
                ScopeImGuiId _id(id++);
                auto& [key, value] = *itr;
                Type<IAutoImGui, _Ty>::GetIType()->DrawAutoImGui(&value, key.c_str(), userdata);
                if (ImGui::Button("erase")) {
                    v.erase(itr);
                    break;
                }
            }
        }
    }
};

template<template<class _Kty, class _Ty, class _Pr, class _Alloc> class ContainerType,
    class _Kty, class _Ty, class _Pr, class _Alloc>
    class Type<IAutoImGui, ContainerType<_Kty, _Ty, _Pr, _Alloc>, std::enable_if_t<std::is_same_v<std::string, _Kty>>>
    : public TypeBase<IAutoImGui, ContainerType<_Kty, _Ty, _Pr, _Alloc>> {
    public:
        using ValueType = ContainerType<_Kty, _Ty, _Pr, _Alloc>;

        struct Userdata : Type< IAutoImGui, _Ty>::Userdata {};

        void DrawAutoImGui(void* addr, const char* name, const UserdataBase* userdata) const override {
            _AutoImGuiMapTypeHelper<ValueType>::DrawAutoImGui(addr, name, userdata);
        }
};

template<template <class _Kty, class _Ty, class _Hasher, class _Keyeq, class _Alloc> class ContainerType,
    class _Kty, class _Ty, class _Hasher, class _Keyeq, class _Alloc>
    class Type<IAutoImGui, ContainerType<_Kty, _Ty, _Hasher, _Keyeq, _Alloc>, std::enable_if_t<std::is_same_v<std::string, _Kty>>>
    : public TypeBase<IAutoImGui, ContainerType<_Kty, _Ty, _Hasher, _Keyeq, _Alloc>> {
    public:
        using ValueType = ContainerType<_Kty, _Ty, _Hasher, _Keyeq, _Alloc>;

        struct Userdata : Type< IAutoImGui, _Ty>::Userdata {};

        void DrawAutoImGui(void* addr, const char* name, const UserdataBase* userdata) const override {
            _AutoImGuiMapTypeHelper<ValueType>::DrawAutoImGui(addr, name, userdata);
        }
};

} // namespace reflection
