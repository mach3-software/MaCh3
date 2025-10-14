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

    #ifdef MaCh3_MEMORY_PROFILING_ENABLED

        // use tracy macros to profile memory allocation and free
        #define MaCh3_ProfileMemory                    \
        void * operator new ( std :: size_t count ) {  \
            auto ptr = malloc ( count ) ;              \
            TracyAlloc ( ptr , count ) ;               \
            return ptr ; }                             \
        void operator delete ( void * ptr ) noexcept { \
            TracyFree ( ptr ) ;                        \
            free ( ptr ) ; }                           \

        
        // use tracy macros to profile memory allocation and free
        #define MaCh3_ProfileMemoryN(name)             \
        void * operator new ( std :: size_t count ) {  \
            auto ptr = malloc ( count ) ;              \
            TracyAllocN ( ptr , count, name ) ;        \
            return ptr ; }                             \
        void operator delete ( void * ptr ) noexcept { \
            TracyFreeN ( ptr, name ) ;                 \
            free ( ptr ) ; }                           \

    #else 

        #define MaCh3_ProfileMemory
        #define MaCh3_ProfileMemoryN

    #endif 

#else
    #define MaCh3_zoneScoped
    #define MaCh3_zoneScopedN(x)  
    #define MaCh3_FrameMark
#endif
