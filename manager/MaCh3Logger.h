#pragma once

#include "spdlog/spdlog.h"

//KS: Based on this https://github.com/gabime/spdlog/blob/a2b4262090fd3f005c2315dcb5be2f0f1774a005/include/spdlog/spdlog.h#L284

#define MACH3LOG_TRACE SPDLOG_TRACE
#define MACH3LOG_DEBUG SPDLOG_DEBUG
#define MACH3LOG_INFO SPDLOG_INFO
#define MACH3LOG_WARN SPDLOG_WARN
#define MACH3LOG_ERROR SPDLOG_ERROR
#define MACH3LOG_CRITICAL SPDLOG_CRITICAL
#define MACH3LOG_OFF SPDLOG_OFF
