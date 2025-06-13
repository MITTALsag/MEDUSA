#pragma once


#ifdef OS_WINDOWS
    // Windows
    #define YES_DLL_EXPORT __declspec(dllexport)
    #define NO_DLL_EXPORT /*no export*/
    //
#elif defined(OS_LINUX)
    // Linux
    #define YES_DLL_EXPORT __attribute__((visibility("default")))
    #define NO_DLL_EXPORT __attribute__((visibility("hidden")))
    //
#else
    // Other OS **HERE**
    // ...
    // ...
    *** Export symbols must be defined! ***
    //
#endif
