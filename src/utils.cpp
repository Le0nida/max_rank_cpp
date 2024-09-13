//
// Created by leona on 13/09/2024.
//

#include "utils.h"

#if defined(_WIN32)
#include <windows.h>

size_t getAvailableMemory() {
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatusEx(&status);
    return static_cast<size_t>(status.ullAvailPhys);
}

#elif defined(__linux__)
#include <sys/sysinfo.h>

size_t getAvailableMemory() {
    struct sysinfo info;
    if (sysinfo(&info) == 0) {
        return static_cast<size_t>(info.freeram) * info.mem_unit;
    }
    return 0;
}

#else
#error "Unsupported platform"
#endif