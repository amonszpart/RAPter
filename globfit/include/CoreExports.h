#pragma once
#ifndef CORE_EXPORTS_H_
#define CORE_EXPORTS_H_

#if defined WIN32 || defined _WIN32 || defined WINCE || defined __MINGW32__
    #ifdef CORE_API_EXPORTS
        #define CORE_EXPORTS __declspec(dllexport)
    #else
        #define CORE_EXPORTS
    #endif
#else
    #define CORE_EXPORTS
#endif

#endif  //#ifndef CORE_EXPORTS_H_
