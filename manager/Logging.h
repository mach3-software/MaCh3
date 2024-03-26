#pragma once

#include <chrono>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <stdarg.h>
#include <ctime>

namespace LOG {
// EM: Define possible log levels
//     Doing it this way allows us to easily add in a new level between or before others
#define MACH3_LOG_kSilent 0      // <- Don't print developer messages, this wont include e.g. progress bars, summary tables etc.
#define MACH3_LOG_kError 1       // <- Print errors
#define MACH3_LOG_kWarning 2     // <- Print warnings
#define MACH3_LOG_kInfo 3        // <- Print high level information
#define MACH3_LOG_kDebug 4       // <- print standard debug information

#define MACH3_LOG_kTrace 5       // <- Print EVERYTHING

// EM: define a preamble to print at the start of log messages
//     This can be used to print the time as is done currently, could also add info about which file or function the message is coming from?
inline void printPreamble(std::string prefix){
    // EM: Print the current time
    const auto now = std::chrono::system_clock::now();

    const std::time_t t_c = std::chrono::system_clock::to_time_t(now);
    const tm *time = localtime(&t_c);
    
    char dateTime[20];
    strftime(dateTime, 100, "%D %T", time);
    
    std::cout << dateTime;

    // EM: Print the prefix
    std::cout << prefix;
}

// fn to log variable number of arguments
// fmt is a format string which uses the same formatting as printf()
// prefix is a prefix to print before the other arguments
// the following arguments are the things to be logged
inline void LOG(std::string prefix, const char*fmt...){
    printPreamble(prefix);
    
    va_list args;
    va_start(args, fmt);

    // EM: print the variable length list of arguments passed
    vprintf(fmt, args);

    std::cout << std::endl;
}

// EM: Define these using macros so they can be stripped from the build depending on the log level so they won't affect performance
#if MACH3_LOG_LEVEL >= MACH3_LOG_kError
    #define ERRORLOG(...)LOG("[\033[1;31mERROR\033[0m]: ", __VA_ARGS__);
#else
    #define ERRORLOG(...)
#endif


#if MACH3_LOG_LEVEL >= MACH3_LOG_kWarning
    #define WARNLOG(...)LOG("[\033[1;33mWARN\033[0m]: ", __VA_ARGS__);
#else
    #define WARNLOG(...)
#endif


#if MACH3_LOG_LEVEL >= MACH3_LOG_kInfo
    #define INFOLOG(...)LOG("[\033[1;92mINFO\033[0m]: ", __VA_ARGS__);
#else
    #define INFOLOG(...)
#endif


#if MACH3_LOG_LEVEL >= MACH3_LOG_kDebug
    #define DEBUGLOG(...)LOG("[\033[1;96mDEBUG\033[0m]: ", __VA_ARGS__);
#else
    #define DEBUGLOG(...)
#endif


#if MACH3_LOG_LEVEL >= MACH3_LOG_kTrace
    #define TRACELOG(...)LOG("[\033[1;90mTRACE\033[0m]: ", __VA_ARGS__);
#else
    #define TRACELOG(...)
#endif

}
