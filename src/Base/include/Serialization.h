#pragma once

#include "log.h"

#define LOG_WARN_MSG(msg) LOG_WARN("{}", msg)
#define FIELD_NOT_FOUND_HANDLE LOG_WARN_MSG
#include <reflection/serialization_ext_glm.h>
#undef LOG_WARN_MSG

using ISerializable = reflection::ISerialization;
using reflection::Serialize;
using reflection::Deserialize;

#define FIELD_DECLARE(name, ...) FIELD_DECLARATION(#name, name, __VA_ARGS__)
#define STRUCT_FIELD_DECLARE(name, ...) STRUCT_FIELD_DECLARATION(#name, name, __VA_ARGS__)
