#pragma once

#include "reflection.h"

#include <memory>
#include <array>
#include <vector>
#include <list>
#include <map>
#include <unordered_map>

#include <rapidjson/document.h>
#include <rapidjson/prettywriter.h>
#include <magic_enum.hpp>

#include "util.h"

#ifndef FIELD_NOT_FOUND_HANDLE
#include <iostream>
#define FIELD_NOT_FOUND_HANDLE(msg) (std::cerr << (msg) << std::endl)
#endif

namespace reflection {

class ISerialization : public IReflectionBase<ISerialization> {
};

template<>
class IType<ISerialization> {
public:
    virtual void Serialize(const void*, rapidjson::PrettyWriter<rapidjson::StringBuffer>&) const = 0;

    virtual void Deserialize(void*, const rapidjson::Value&) const = 0;
};

template<class T>
void Serialize(const T& object, rapidjson::PrettyWriter<rapidjson::StringBuffer>& writer) {
    Type<ISerialization, T>::GetIType()->Serialize(&object, writer);
}

template<class T>
void Deserialize(T& object, const rapidjson::Value& value) {
    Type<ISerialization, T>::GetIType()->Deserialize(&object, value);
}

template<>
class Type<ISerialization, int> : public TypeBase<ISerialization, int> {
public:
    void Serialize(const void* addr, rapidjson::PrettyWriter<rapidjson::StringBuffer>& writer) const override {
        const auto& v = *static_cast<const ValueType*>(addr);
        writer.Int(v);
    }

    void Deserialize(void* addr, const rapidjson::Value& value) const override {
        R_ASSERT(value.IsInt());
        auto& v = *static_cast<ValueType*>(addr);
        v = value.GetInt();
    }
};

template<>
class Type<ISerialization, bool> : public TypeBase<ISerialization, bool> {
public:
    void Serialize(const void* addr, rapidjson::PrettyWriter<rapidjson::StringBuffer>& writer) const override {
        const auto& v = *static_cast<const ValueType*>(addr);
        writer.Bool(v);
    }

    void Deserialize(void* addr, const rapidjson::Value& value) const override {
        R_ASSERT(value.IsBool());
        auto& v = *static_cast<ValueType*>(addr);
        v = value.GetBool();
    }
};

template<>
class Type<ISerialization, float> : public TypeBase<ISerialization, float> {
public:
    void Serialize(const void* addr, rapidjson::PrettyWriter<rapidjson::StringBuffer>& writer) const override {
        const auto& v = *static_cast<const ValueType*>(addr);
        writer.Double(static_cast<double>(v));
    }

    void Deserialize(void* addr, const rapidjson::Value& value) const override {
        R_ASSERT(value.IsNumber());
        auto& v = *static_cast<ValueType*>(addr);
        v = value.GetFloat();
    }
};

template<>
class Type<ISerialization, double> : public TypeBase<ISerialization, double> {
public:
    void Serialize(const void* addr, rapidjson::PrettyWriter<rapidjson::StringBuffer>& writer) const override {
        const auto& v = *static_cast<const ValueType*>(addr);
        writer.Double(v);
    }

    void Deserialize(void* addr, const rapidjson::Value& value) const override {
        R_ASSERT(value.IsNumber());
        auto& v = *static_cast<ValueType*>(addr);
        v = value.GetDouble();
    }
};

template<>
class Type<ISerialization, std::string> : public TypeBase<ISerialization, std::string> {
public:
    void Serialize(const void* addr, rapidjson::PrettyWriter<rapidjson::StringBuffer>& writer) const override {
        const auto& v = *static_cast<const ValueType*>(addr);
        writer.String(v.c_str());
    }

    void Deserialize(void* addr, const rapidjson::Value& value) const override {
        R_ASSERT(value.IsString());
        auto& v = *static_cast<ValueType*>(addr);
        v = value.GetString();
    }
};

template<class T>
class Type<ISerialization, T, std::enable_if_t<std::is_enum_v<T>>> : public TypeBase<ISerialization, T> {
public:
    static_assert(!std::is_same_v<ISerialization, T>);

    void Serialize(const void* addr, rapidjson::PrettyWriter<rapidjson::StringBuffer>& writer) const override {
        const auto v = *static_cast<const T*>(addr);
        writer.String(std::string(magic_enum::enum_name(v)).c_str());
    }

    void Deserialize(void* addr, const rapidjson::Value& value) const override {
        R_ASSERT(value.IsString());
        auto e = magic_enum::enum_cast<T>(value.GetString());
        R_ASSERT(e.has_value());
        auto& v = *static_cast<T*>(addr);
        v = e.value();
    }
};

template<class T>
class Type<ISerialization, T, std::enable_if_t<std::is_base_of_v<ISerialization, T> || IsReflectableStruct<ISerialization, T>::value>>
    : public TypeBase<ISerialization, T> {
public:
    using ValueType = T;

    void Serialize(const void* addr, rapidjson::PrettyWriter<rapidjson::StringBuffer>& writer) const override {
        const auto& v = *static_cast<const ValueType*>(addr);
        writer.StartObject();
        for (const auto& [name, fun] : GetFieldTable(v, static_cast<ISerialization*>(nullptr))) {
            writer.String(name.c_str());
            auto info = fun(const_cast<ValueType*>(&v));
            info.type->Serialize(info.address, writer);
        }
        writer.EndObject();
    }

    void Deserialize(void* addr, const rapidjson::Value& value) const override {
        R_ASSERT(value.IsObject());
        auto& v = *static_cast<ValueType*>(addr);
        for (const auto& [name, fun] : GetFieldTable(v, static_cast<ISerialization*>(nullptr))) {
            auto itr = value.FindMember(name.c_str());
            if (itr != value.MemberEnd()) {
                auto info = fun(&v);
                info.type->Deserialize(info.address, itr->value);
            }
            else {
                FIELD_NOT_FOUND_HANDLE("Field \"" + name + "\" not found");
            }
        }
    }
};

template<class _Ty, class _Dx>
class Type<ISerialization, std::unique_ptr<_Ty, _Dx>, std::enable_if_t<std::is_base_of_v<ISerialization, _Ty>&& SubclassInfo<_Ty>::has>>
    : public TypeBase<ISerialization, std::unique_ptr<_Ty, _Dx>> {
public:
    using ValueType = std::unique_ptr<_Ty, _Dx>;
    static constexpr auto kTypeKey = "type";
    static constexpr auto kDataKey = "data";

    void Serialize(const void* addr, rapidjson::PrettyWriter<rapidjson::StringBuffer>& writer) const override {
        const auto& v = *static_cast<const ValueType*>(addr);
        if (v) {
            writer.StartObject();
            writer.String(kTypeKey);
            writer.String(typeid(*v).name());
            writer.String(kDataKey);
            Type<ISerialization, _Ty>::GetIType()->Serialize(v.get(), writer);
            writer.EndObject();
        }
        else {
            writer.Null();
        }
    }

    void Deserialize(void* addr, const rapidjson::Value& value) const override {
        R_ASSERT(value.IsObject() || value.IsNull());
        auto& v = *static_cast<ValueType*>(addr);
        if (value.IsObject()) {
            auto typeitr = value.FindMember(kTypeKey);
            R_ASSERT(typeitr != value.MemberEnd());
            auto& name = typeitr->value;
            R_ASSERT(name.IsString());
            v.reset(SubclassInfo<_Ty>::GetFactoryTable().at(name.GetString())());

            auto dataitr = value.FindMember(kDataKey);
            R_ASSERT(dataitr != value.MemberEnd());
            Type<ISerialization, _Ty>::GetIType()->Deserialize(v.get(), dataitr->value);
        }
    }
};

template<class _Ty, class _Dx>
class Type<ISerialization, std::unique_ptr<_Ty, _Dx>, std::enable_if_t<!(std::is_base_of_v<ISerialization, _Ty>&& SubclassInfo<_Ty>::has)>>
    : public TypeBase<ISerialization, std::unique_ptr<_Ty, _Dx>> {
public:
    using ValueType = std::unique_ptr<_Ty, _Dx>;

    void Serialize(const void* addr, rapidjson::PrettyWriter<rapidjson::StringBuffer>& writer) const override {
        const auto& v = *static_cast<const ValueType*>(addr);
        if (v)
            Type<ISerialization, _Ty>::GetIType()->Serialize(v.get(), writer);
        else
            writer.Null();
    }

    void Deserialize(void* addr, const rapidjson::Value& value) const override {
        auto& v = *static_cast<ValueType*>(addr);
        if (value.IsNull()) {
            v.reset();
        }
        else {
            if (v == nullptr)
                v.reset(new _Ty());
            Type<ISerialization, _Ty>::GetIType()->Deserialize(v.get(), value);
        }
    }
};

template<class _Ty, size_t _Size>
struct _SerializationArrayTypeHelper {
    static void Serialize(const void* addr, rapidjson::PrettyWriter<rapidjson::StringBuffer>& writer) {
        const auto arr = static_cast<const _Ty*>(addr);
        writer.StartArray();
        for (size_t i = 0; i < _Size; ++i) {
            Type<ISerialization, _Ty>::GetIType()->Serialize(&arr[i], writer);
        }
        writer.EndArray();
    }

    static void Deserialize(void* addr, const rapidjson::Value& value) {
        R_ASSERT(value.IsArray());
        R_ASSERT(value.Size() == _Size);
        auto arr = static_cast<_Ty*>(addr);
        for (size_t i = 0; i < _Size; ++i) {
            Type<ISerialization, _Ty>::GetIType()->Deserialize(&arr[i], value[static_cast<rapidjson::SizeType>(i)]);
        }
    }
};

template<class _Ty, size_t _Size>
class Type< ISerialization, _Ty[_Size]>
    : public TypeBase<ISerialization, _Ty[_Size]> {
public:
    void Serialize(const void* addr, rapidjson::PrettyWriter<rapidjson::StringBuffer>& writer) const override {
        _SerializationArrayTypeHelper<_Ty, _Size>::Serialize(addr, writer);
    }

    void Deserialize(void* addr, const rapidjson::Value& value) const override {
        _SerializationArrayTypeHelper<_Ty, _Size>::Deserialize(addr, value);
    }
};

template<class _Ty, size_t _Size>
class Type< ISerialization, std::array<_Ty, _Size>>
    : public TypeBase<ISerialization, std::array<_Ty, _Size>> {
public:
    void Serialize(const void* addr, rapidjson::PrettyWriter<rapidjson::StringBuffer>& writer) const override {
        _SerializationArrayTypeHelper<_Ty, _Size>::Serialize(addr, writer);
    }

    void Deserialize(void* addr, const rapidjson::Value& value) const override {
        _SerializationArrayTypeHelper<_Ty, _Size>::Deserialize(addr, value);
    }
};

template<template<class _Ty, class _Alloc> class ContainerType, class _Ty, class _Alloc>
class Type<ISerialization, ContainerType<_Ty, _Alloc>
    , std::enable_if_t<std::is_same_v<ContainerType<_Ty, _Alloc>, std::vector<_Ty, _Alloc>>
                    || std::is_same_v<ContainerType<_Ty, _Alloc>, std::list<_Ty, _Alloc>>>>
    : public TypeBase<ISerialization, ContainerType<_Ty, _Alloc>> {
public:
    using ValueType = ContainerType<_Ty, _Alloc>;

    void Serialize(const void* addr, rapidjson::PrettyWriter<rapidjson::StringBuffer>& writer) const override {
        const auto& v = *static_cast<const ValueType*>(addr);
        writer.StartArray();
        for (const auto& e : v) {
            Type<ISerialization, _Ty>::GetIType()->Serialize(&e, writer);
        }
        writer.EndArray();
    }

    void Deserialize(void* addr, const rapidjson::Value& value) const override {
        R_ASSERT(value.IsArray());
        auto& v = *static_cast<ValueType*>(addr);

        v.clear();
        for (const auto& e : value.GetArray()) {
            _Ty tmp{};
            Type<ISerialization, _Ty>::GetIType()->Deserialize(&tmp, e);
            v.emplace_back(std::move(tmp));
        }
    }
};

template<class T>
struct _SerializationMapTypeHelper {
    using _Ty = typename T::mapped_type;

    static void Serialize(const void* addr, rapidjson::PrettyWriter<rapidjson::StringBuffer>& writer) {
        const auto& v = *static_cast<const T*>(addr);
        writer.StartObject();
        for (const auto& [key, value] : v) {
            writer.String(key.c_str());
            Type<ISerialization, _Ty>::GetIType()->Serialize(&value, writer);
        }
        writer.EndObject();
    }

    static void Deserialize(void* addr, const rapidjson::Value& value) {
        R_ASSERT(value.IsObject());
        auto& v = *static_cast<T*>(addr);
        
        v.clear();
        for (const auto& e : value.GetObject()) {
            _Ty tmp{};
            Type<ISerialization, _Ty>::GetIType()->Deserialize(&tmp, e.value);
            v.emplace(e.name.GetString(), std::move(tmp));
        }
    }
};

template<template<class _Kty, class _Ty, class _Pr, class _Alloc> class ContainerType,
    class _Kty, class _Ty, class _Pr, class _Alloc>
    class Type<ISerialization, ContainerType<_Kty, _Ty, _Pr, _Alloc>, std::enable_if_t<std::is_same_v<std::string, _Kty>>>
    : public TypeBase<ISerialization, ContainerType<_Kty, _Ty, _Pr, _Alloc>> {
    public:
        using ValueType = ContainerType<_Kty, _Ty, _Pr, _Alloc>;

        void Serialize(const void* addr, rapidjson::PrettyWriter<rapidjson::StringBuffer>& writer) const override {
            _SerializationMapTypeHelper<ValueType>::Serialize(addr, writer);
        }

        void Deserialize(void* addr, const rapidjson::Value& value) const override {
            _SerializationMapTypeHelper<ValueType>::Deserialize(addr, value);
        }
};

template<template <class _Kty, class _Ty, class _Hasher, class _Keyeq, class _Alloc> class ContainerType,
    class _Kty, class _Ty, class _Hasher, class _Keyeq, class _Alloc>
    class Type<ISerialization, ContainerType<_Kty, _Ty, _Hasher, _Keyeq, _Alloc>, std::enable_if_t<std::is_same_v<std::string, _Kty>>>
    : public TypeBase<ISerialization, ContainerType<_Kty, _Ty, _Hasher, _Keyeq, _Alloc>> {
    public:
        using ValueType = ContainerType<_Kty, _Ty, _Hasher, _Keyeq, _Alloc>;

        void Serialize(const void* addr, rapidjson::PrettyWriter<rapidjson::StringBuffer>& writer) const override {
            _SerializationMapTypeHelper<ValueType>::Serialize(addr, writer);
        }

        void Deserialize(void* addr, const rapidjson::Value& value) const override {
            _SerializationMapTypeHelper<ValueType>::Deserialize(addr, value);
        }
};

} // namespace reflection
