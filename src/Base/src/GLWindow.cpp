#include "GLWindow.h"

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <memory>

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <stb_image_write.h>

#include "Utils.h"
#include "gl.hpp"
#include "log.h"
#include "Profiler.h"

GLWindow::GLWindow(const char* name, int width, int height, bool vsync) {
    SetCurrentDirToExe(); // 切换到exe所在目录

    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#if _DEBUG
    glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GLFW_TRUE);
#endif

    window = glfwCreateWindow(width, height, name, NULL, NULL);
    if (window == NULL)
    {
        throw std::runtime_error("Failed to create GLFW window");
        glfwTerminate();
    }
    glfwMakeContextCurrent(window);
    if (!vsync)
        glfwSwapInterval(0);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        throw std::runtime_error("Failed to initialize GLAD");
    }

#if _DEBUG
    // https://learnopengl.com/In-Practice/Debugging
    int flags; glGetIntegerv(GL_CONTEXT_FLAGS, &flags);
    if (flags & GL_CONTEXT_FLAG_DEBUG_BIT)
    {
        glEnable(GL_DEBUG_OUTPUT);
        glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS); // makes sure errors are displayed synchronously
        glDebugMessageCallback([](GLenum source,
            GLenum type,
            unsigned int id,
            GLenum severity,
            GLsizei length,
            const char* message,
            const void* userParam)
            {
                // ignore non-significant error/warning codes
                if (id == 131169 || id == 131185 || id == 131218 || id == 131204) return;

                std::stringstream msg;
                msg << "---------------\n";
                msg << "Debug message (0x" << std::hex << id << "): \t";

                switch (type)
                {
                case GL_DEBUG_TYPE_ERROR:
                    msg << "Type: Error"; break;
                case GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR:
                    msg << "Type: Deprecated Behaviour"; break;
                case GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR:
                    msg << "Type: Undefined Behaviour"; break;
                case GL_DEBUG_TYPE_PORTABILITY:
                    msg << "Type: Portability"; break;
                case GL_DEBUG_TYPE_PERFORMANCE:
                    msg << "Type: Performance"; break;
                case GL_DEBUG_TYPE_MARKER:
                    msg << "Type: Marker"; break;
                case GL_DEBUG_TYPE_PUSH_GROUP:
                    msg << "Type: Push Group"; return; // Do not print
                case GL_DEBUG_TYPE_POP_GROUP:
                    msg << "Type: Pop Group"; return; // Do not print
                case GL_DEBUG_TYPE_OTHER:
                    msg << "Type: Other"; break;
                }
                msg << "\t";

                switch (source)
                {
                case GL_DEBUG_SOURCE_API:
                    msg << "Source: API"; break;
                case GL_DEBUG_SOURCE_WINDOW_SYSTEM:
                    msg << "Source: Window System"; break;
                case GL_DEBUG_SOURCE_SHADER_COMPILER:
                    msg << "Source: Shader Compiler"; break;
                case GL_DEBUG_SOURCE_THIRD_PARTY:
                    msg << "Source: Third Party"; break;
                case GL_DEBUG_SOURCE_APPLICATION:
                    msg << "Source: Application"; break;
                case GL_DEBUG_SOURCE_OTHER:
                    msg << "Source: Other"; break;
                }
                msg << "\t";

                switch (severity)
                {
                case GL_DEBUG_SEVERITY_HIGH:
                    msg << "Severity: high"; break;
                case GL_DEBUG_SEVERITY_MEDIUM:
                    msg << "Severity: medium"; break;
                case GL_DEBUG_SEVERITY_LOW:
                    msg << "Severity: low"; break;
                case GL_DEBUG_SEVERITY_NOTIFICATION:
                    msg << "Severity: notification"; break;
                }
                msg << "\n";
                msg << message;
                LOG_WARN("{}", msg.str());
            }, nullptr);
        glDebugMessageControl(GL_DONT_CARE, GL_DONT_CARE, GL_DONT_CARE, 0, nullptr, GL_TRUE);
    }
    else {
        LOG_ERROR("Failed to create OpenGL Debug Context");
    }
#endif

    glfwSetWindowUserPointer(window, this);
    glfwSetFramebufferSizeCallback(window,
        [](GLFWwindow* window, int width, int height) {
            auto w = static_cast<GLWindow*>(glfwGetWindowUserPointer(window));
            w->HandleReshapeEvent(width, height);
        });
    glfwSetKeyCallback(window,
        [](GLFWwindow* window, int key, int scancode, int action, int mods) {
            if (ImGui::GetIO().WantCaptureKeyboard)
                return;
            auto w = static_cast<GLWindow*>(glfwGetWindowUserPointer(window));
            if (action == GLFW_PRESS)
                w->HandleKeyboardEvent(key);
        });
    glfwSetCursorPosCallback(window,
        [](GLFWwindow* window, double x, double y) {
            if (ImGui::GetIO().WantCaptureMouse)
                return;
            auto w = static_cast<GLWindow*>(glfwGetWindowUserPointer(window));
            w->HandleMouseEvent(x, y);
        });
    glfwSetScrollCallback(window,
        [](GLFWwindow* window, double xoffset, double yoffset) {
            auto w = static_cast<GLWindow*>(glfwGetWindowUserPointer(window));
            w->HandleMouseWheelEvent(yoffset);
        });
    glfwSetDropCallback(window,
        [](GLFWwindow* window, int count, const char** paths) {
            auto w = static_cast<GLWindow*>(glfwGetWindowUserPointer(window));
            w->HandleDropEvent(count, paths);
        });

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;
    io.ConfigFlags |= ImGuiConfigFlags_ViewportsEnable;
    ImGui::StyleColorsDark();
    ImGuiStyle& style = ImGui::GetStyle();
    if (io.ConfigFlags & ImGuiConfigFlags_ViewportsEnable)
    {
        style.WindowRounding = 0.0f;
        style.Colors[ImGuiCol_WindowBg].w = 1.0f;
    }
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 460");
}

GLWindow::~GLWindow() {
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwTerminate();
}

void GLWindow::MainLoop() {
    while (!glfwWindowShouldClose(window)) {
        Profiler::Instance().NewFrame();
        {
            PERF_MARKER("Frame")
            HandleDisplayEvent();
            {
                PERF_MARKER("GUI")
                ImGui_ImplOpenGL3_NewFrame();
                ImGui_ImplGlfw_NewFrame();
                ImGui::NewFrame();

                HandleDrawGuiEvent();
                CheckError();

                ImGui::Render();
                ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
            }
        }
        if (ImGui::GetIO().ConfigFlags & ImGuiConfigFlags_ViewportsEnable) {
            GLFWwindow* backup_current_context = glfwGetCurrentContext();
            ImGui::UpdatePlatformWindows();
            ImGui::RenderPlatformWindowsDefault();
            glfwMakeContextCurrent(backup_current_context);
        }
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
}

std::tuple<int, int>  GLWindow::GetWindowSize() {
    int width, height;
    glfwGetWindowSize(window, &width, &height);
    return { width, height };
}

void GLWindow::ResizeWindow(int width, int height) {
    glfwSetWindowSize(window, width, height);
}

void GLWindow::SetFullScreen(bool is_full_screen) {
    if (is_full_screen) {
        glfwGetWindowPos(window, &xpos_, &ypos_);
        glfwGetWindowSize(window, &width_, &height_);
        const GLFWvidmode* mode = glfwGetVideoMode(glfwGetPrimaryMonitor());
        glfwSetWindowMonitor(window, glfwGetPrimaryMonitor(), 0, 0, mode->width, mode->height, 0);
    }
    else {
        // restore last window size and position
        glfwSetWindowMonitor(window, nullptr, xpos_, ypos_, width_, height_, 0);
    }
}

void GLWindow::Close() {
    glfwSetWindowShouldClose(window, true);
}

void GLWindow::Error(const std::string& msg) {
    b_current_frame_error_ = true;
    if (err_msgs_.size() >= 60)
        err_msgs_.pop_front();
    err_msgs_.push_back(msg);
    LOG_ERROR("{}", msg);
}

void GLWindow::ScreenShot(const std::string& path_no_ext, ImageFormat format) {
    auto [width, height] = GetWindowSize();
    auto pixels = std::make_unique<std::byte[]>(static_cast<size_t>(width) * height * 3);
    glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels.get());
    switch (format) {
    case GLWindow::PNG:
        stbi_write_png((path_no_ext + ".png").c_str(), width, height, 3, pixels.get(), width * 3);
        break;
    case GLWindow::JPG:
        stbi_write_jpg((path_no_ext + ".jpg").c_str(), width, height, 3, pixels.get(), 90);
        break;
    default:
        break;
    }
}

void GLWindow::CheckError() {
    if (b_current_frame_error_)
        ImGui::OpenPopup("Error");
    if (ImGui::BeginPopupModal("Error", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) {
        b_current_frame_error_ = false;
        for (const auto& msg : err_msgs_)
            ImGui::Text(msg.c_str());
        if (ImGui::Button("OK", ImVec2(120, 0))) {
            err_msgs_.clear();
            ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
    }
}
