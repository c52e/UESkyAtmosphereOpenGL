#pragma once

#include <array>

#include <glad/glad.h>

template<class T, GLsizei N>
class UniqueHandles {
public:
    UniqueHandles() {
        ids_.fill(0);
    }

    ~UniqueHandles() {
        T::Delete(N, ids_.data());
    }

    UniqueHandles(const UniqueHandles&) = delete;

    UniqueHandles(UniqueHandles&& rhs) noexcept
        : UniqueHandles() {
        swap(rhs);
    }

    UniqueHandles& operator=(UniqueHandles rhs) noexcept {
        swap(rhs);
        return *this;
    }

    template<class = std::enable_if_t<!T::create_with_param>>
    void Create() {
        T::Create(N, ids_.data());
    }

    template<class = std::enable_if_t<T::create_with_param>>
    void Create(GLenum target) {
        T::Create(target, N, ids_.data());
    }

    const GLuint& operator[](GLsizei i) const noexcept {
        return ids_[i];
    }

    constexpr auto cbegin() const noexcept {
        return ids_.cbegin();
    }

    constexpr auto cend() const noexcept {
        return ids_.cend();
    }

    constexpr auto crbegin() const noexcept {
        return ids_.crbegin();
    }

    constexpr auto crend() const noexcept {
        return ids_.crend();
    }

    constexpr GLsizei size() const noexcept {
        return static_cast<GLsizei>(ids_.size());
    }

    void swap(UniqueHandles& rhs) noexcept {
        ids_.swap(rhs.ids_);
    }
private:
    std::array<GLuint, N> ids_;
};

template<class T>
class UniqueHandle {
public:
    template<class = std::enable_if_t<!T::create_with_param>>
    void Create() {
        impl_.Create();
    }

    template<class = std::enable_if_t<T::create_with_param>>
    void Create(GLenum target) {
        impl_.Create(target);
    }

    GLuint id() const noexcept {
        return impl_[0];
    }
private:
    UniqueHandles<T, 1> impl_;
};

#define DECLARE_USING(name, names) \
using GL##name = UniqueHandle<name##Traits>; \
template<GLsizei N> \
using GL##names = UniqueHandles<name##Traits, N>;

#define DECLARE_TRAITS(name, names) \
struct name##Traits { \
static constexpr bool create_with_param = false; \
static void Create(GLsizei n, GLuint* ids) { glCreate##names(n, ids); } \
static void Delete(GLsizei n, GLuint* ids) { glDelete##names(n, ids); } }; \
DECLARE_USING(name, names)

DECLARE_TRAITS(VertexArray, VertexArrays)
DECLARE_TRAITS(Buffer, Buffers)
DECLARE_TRAITS(Framebuffer, Framebuffers)
DECLARE_TRAITS(Sampler, Samplers)

#undef DECLARE_TRAITS

#define DECLARE_TRAITS(name, names) \
struct name##Traits { \
static constexpr bool create_with_param = true; \
static void Create(GLenum target, GLsizei n, GLuint* ids) { glCreate##names(target, n, ids); } \
static void Delete(GLsizei n, GLuint* ids) { glDelete##names(n, ids); } }; \
DECLARE_USING(name, names)

DECLARE_TRAITS(Texture, Textures)
DECLARE_TRAITS(Query, Queries)

#undef DECLARE_TRAITS
#undef DECLARE_USING


#define DEFINE_BIND_UNITS(name) \
template<GLsizei N> inline void GLBind##name(const GLuint (&arr)[N], GLuint first=0) { glBind##name(first, N, arr); }

DEFINE_BIND_UNITS(Textures)
DEFINE_BIND_UNITS(Samplers)
DEFINE_BIND_UNITS(ImageTextures)

#undef DEFINE_BIND_UNITS


#define CONCAT(a, b) CONCAT_INNER(a, b)
#define CONCAT_INNER(a, b) a ## b

#define PERF_MARKER(message) DebugGroup CONCAT(__debug_group__, __LINE__)(message);

class DebugGroup {
public:
    DebugGroup(const char* message) {
        glPushDebugGroup(GL_DEBUG_SOURCE_APPLICATION, 0, -1, message);
    }
    ~DebugGroup() {
        glPopDebugGroup();
    }
};
