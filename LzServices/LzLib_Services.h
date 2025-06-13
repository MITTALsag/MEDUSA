#pragma once

#include <LzLib/LzLib.h>

#ifdef EXPORT_LzServices
    #define DLL_EXPORT_LzServices YES_DLL_EXPORT
#else
    #define DLL_EXPORT_LzServices NO_DLL_EXPORT
#endif
