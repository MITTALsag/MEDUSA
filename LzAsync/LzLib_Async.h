#pragma once

#include <LzLib/LzLib.h>

#ifdef EXPORT_LzAsync
    #define DLL_EXPORT_LzAsync YES_DLL_EXPORT
#else
    #define DLL_EXPORT_LzAsync NO_DLL_EXPORT
#endif
