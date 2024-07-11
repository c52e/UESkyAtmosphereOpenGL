#include "log.h"

static const char* kLogLevelStrings[] = {
	"debug",
	"info",
	"warn",
	"error",
	"fatal",
};

void log(const LogLevel level, const std::string_view message
	, const std::source_location source) {
	const auto levelstr = kLogLevelStrings[static_cast<uint32_t>(level)];
	std::cout << std::format("[{}] {}", levelstr, message) << std::endl;
}

void assert_inner(const std::string_view expression
	, const std::source_location source) {
	LOG_FATAL("Assertion failed: ({}) at file {} (line {})", 
		expression, source.file_name(), source.line());
	exit(0);
}
