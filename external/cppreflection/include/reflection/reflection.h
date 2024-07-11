#pragma once

#include <type_traits>
#include <map>
#include <string>

// Declare in class definition
#define FIELD_DECLARATION_BEGIN(interfacename)                                    \
    const reflection::IReflectionBase<interfacename>::FieldTable& GetFieldTable(  \
            interfacename* = nullptr) const override {                            \
        using Class = std::decay_t<decltype(*this)>;                              \
        using InterfaceType = interfacename;                                      \
        static reflection::IReflectionBase<InterfaceType>::FieldTable m{

#define FIELD_DECLARATION(name, field, ...)                                       \
    { name, [](InterfaceType* arg) {                                              \
        auto p = static_cast<Class*>(arg);                                        \
        using T = decltype(p->field);                                             \
        static reflection::Type<InterfaceType, T>::Userdata d;                    \
        __VA_ARGS__;                                                              \
        return reflection::IReflectionBase<InterfaceType>::FieldInfo {            \
            static_cast<void*>(&(p->field)),                                      \
            reflection::Type<InterfaceType, T>::GetIType(),                       \
            static_cast<reflection::UserdataBase*>(&d) }; }},

#define FIELD_DECLARATION_END_WITH_BASE_CLASS(BaseClass) };  \
        static bool initialized = false;                     \
        if (!initialized) {                                  \
            const auto& base = BaseClass::GetFieldTable(     \
                static_cast<InterfaceType*>(nullptr));       \
            m.insert(base.begin(), base.end());              \
            initialized = true;}                             \
        return m; }


#define FIELD_DECLARATION_END() }; return m; }

// Declare in base class header file (after class definition)
#define HAS_SUBCLASS(classname)                                  \
    template<>                                                   \
    struct reflection::SubclassInfo<classname> {                 \
        using Class = classname;                                 \
        using FactoryFunc = Class * (*) ();                      \
        using FactoryTable = std::map<std::string, FactoryFunc>; \
        static constexpr bool has = true;                        \
        static const FactoryTable& GetFactoryTable();            \
    };

// Declare in base class cpp file
// Subclass header should be included
#define SUBCLASS_DECLARATION_BEGIN(classname)                    \
    const reflection::SubclassInfo<classname>::FactoryTable&     \
        reflection::SubclassInfo<classname>::GetFactoryTable() { \
            static const FactoryTable m {

#define SUBCLASS_DECLARATION(subclass)                   \
    { typeid(subclass).name(), []()                      \
        { return static_cast<Class*>(new subclass()); }},

#define SUBCLASS_DECLARATION_END()  }; return m; }


// Declare in struct header file (after struct definition)
#define STRUCT_FIELD_DECLARATION_BEGIN(structname, interfacename)                  \
    template<> struct reflection::IsReflectableStruct<interfacename, structname> { \
        static constexpr bool value = true;  };                                    \
    template<> inline const auto& reflection::GetFieldTable(                       \
            const structname& object, interfacename* ) {                           \
        using Class = structname;                                                  \
        using InterfaceType = interfacename;                                       \
        using FieldInfo = reflection::IReflectionBase<InterfaceType>::FieldInfo;   \
        using GetFieldFunc = FieldInfo(*) (Class*);                                \
        static std::map<std::string, GetFieldFunc> m{

#define STRUCT_FIELD_DECLARATION(name, field, ...)                    \
    { name, [](Class* p) {                                            \
        using T = decltype(p->field);                                 \
        static reflection::Type<InterfaceType, T>::Userdata d;        \
        __VA_ARGS__;                                                  \
        return FieldInfo {                                            \
            static_cast<void*>(&(p->field)),                          \
            reflection::Type<InterfaceType, T>::GetIType(),           \
            static_cast<reflection::UserdataBase*>(&d) }; }},

#define STRUCT_FIELD_DECLARATION_END()  }; return m; }

namespace reflection {

template<class I>
class IType;

struct UserdataBase {};

template <class I>
class IReflectionBase {
public:
    struct FieldInfo {
        void* address;
        const IType<I>* type;
        const UserdataBase* userdata;
    };

    using GetFieldFunc = FieldInfo(*) (I*);
    using FieldTable = std::map<std::string, GetFieldFunc>;

    // Arg is used to avoid signature collision when a class implements multiple interfaces
    virtual const FieldTable& GetFieldTable(I* = nullptr) const = 0;

    virtual ~IReflectionBase() {};
};

template<class I, class T>
const auto& GetFieldTable(const T& object, I* = nullptr) {
    return object.GetFieldTable(static_cast<I*>(nullptr));
}

template<class I, class T, class Enable = void>
class Type;

template<class I, class T>
class TypeBase : public IType<I> {
public:
    using ValueType = T;

    struct Userdata : UserdataBase {};

    static const IType<I>* GetIType() {
        static const Type<I, T> instance;
        return &instance;
    }
protected:
    TypeBase() {};
};


template<class T>
struct SubclassInfo {
    static constexpr bool has = false;
};

template<class I, class T>
struct IsReflectableStruct {
    static constexpr bool value = false;
};

} // namespace reflection
