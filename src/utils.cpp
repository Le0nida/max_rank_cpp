//
// Created by leona on 13/09/2024.
//

#include "utils.h"

#include <iostream>

#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>

size_t getAvailableMemory() {
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatusEx(&status);
    return static_cast<size_t>(status.ullAvailPhys);
}

void printMemoryUsage(const std::string& descr) {
    PROCESS_MEMORY_COUNTERS memInfo;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &memInfo, sizeof(memInfo))) {
        const SIZE_T virtualMemUsedByMe = memInfo.WorkingSetSize;
        if (!descr.empty())
        {
            std::cout << descr << " - ";
        }
        std::cout << "Memory Usage (WorkingSetSize): " << virtualMemUsedByMe / (1024 * 1024) << " MB" << std::endl;
    } else {
        std::cerr << "Failed to get memory info" << std::endl;
    }
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