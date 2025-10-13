#pragma once 

#include <iostream>

#ifdef TRACY_ENABLE
    #include <tracy/Tracy.hpp>

    #define TRACY_SAMPLING_HZ 1000e3

    #ifdef MaCh3_INSTRUMENTED_PROFILING_ENABLED
        #define MaCh3_zoneScoped ZoneScoped 
        #define MaCh3_zoneScopedN(x) ZoneScopedN(x) 

    #else 
        #define MaCh3_zoneScoped TracyNoop
        #define MaCh3_zoneScopedN(x) ZoneScopedN

    #endif

    #define MaCh3_FrameMark FrameMark

#else
    #define MaCh3_zoneScoped
    #define MaCh3_zoneScopedN(x)  
    #define MaCh3_FrameMark
#endif
