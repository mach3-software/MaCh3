#pragma once

#include <iostream>

#ifdef TRACY_ENABLE

    #define TRACY_DELAYED_INIT

    #include <tracy/Tracy.hpp>

    #define TRACY_SAMPLING_HZ 1000e3

    #ifdef MaCh3_INSTRUMENTED_PROFILING_ENABLED
        #define MaCh3_ProfileScope ZoneScoped 
        #define MaCh3_ProfileScopeN(x) ZoneScopedN(x) 

    #else 
        #define MaCh3_ProfileScope TracyNoop
        #define MaCh3_ProfileScopeN(x) ZoneScopedN

    #endif

    #define MaCh3_FrameMark FrameMark

    #ifdef MaCh3_MEMORY_PROFILING_ENABLED

        // use tracy macros to profile memory allocation and free
        #define MaCh3_ProfileMemory                    \
        inline void * operator new ( std :: size_t count ) {  \
            auto ptr = malloc ( count ) ;              \
            TracyAlloc ( ptr , count ) ;               \
            return ptr ; }                             \
        inline void operator delete ( void * ptr ) noexcept { \
            TracyFree ( ptr ) ;                        \
            free ( ptr ) ; }                           \

        
        // use tracy macros to profile memory allocation and free
        #define MaCh3_ProfileMemoryN(name)             \
        inline void * operator new ( std :: size_t count ) {  \
            auto ptr = malloc ( count ) ;              \
            TracyAllocN ( ptr , count, name ) ;        \
            return ptr ; }                             \
        inline void operator delete ( void * ptr ) noexcept { \
            TracyFreeN ( ptr, name ) ;                 \
            free ( ptr ) ; }                           \

    #else 

        #define MaCh3_ProfileMemory
        #define MaCh3_ProfileMemoryN(name)

    #endif 

#else
    #define MaCh3_ProfileScope
    #define MaCh3_ProfileScopeN(x)  
    #define MaCh3_FrameMark
    #define MaCh3_ProfileMemory
    #define MaCh3_ProfileMemoryN(name)
#endif
