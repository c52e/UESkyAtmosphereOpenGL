#pragma once

#include <iostream>
#include <format>
#include <source_location>

enum class LogLevel {
    Debug = 0,
    Info,
    Warn,
    Error,
    Fatal,
};

void log(const LogLevel level, const std::string_view message
    , const std::source_location source = std::source_location::current());

void assert_inner(const std::string_view expression
    , const std::source_location source = std::source_location::current());

#define LOG_DEBUG(...) log(LogLevel::Debug, std::format(__VA_ARGS__));
#define LOG_INFO(...) log(LogLevel::Info, std::format(__VA_ARGS__));
#define LOG_WARN(...) log(LogLevel::Warn, std::format(__VA_ARGS__));
#define LOG_ERROR(...) log(LogLevel::Error, std::format(__VA_ARGS__));
#define LOG_FATAL(...) log(LogLevel::Fatal, std::format(__VA_ARGS__));

#define ASSERT(expression) (void)((!!(expression)) || (assert_inner(#expression), 0))
