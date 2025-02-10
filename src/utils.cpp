#include "utils.h"
#include <cstdio>
#include <cstdlib>   // for exit(0)
#include <cstring>   // for strtok
#include <cmath>     // for std::ceil
#include <filesystem>

#if defined(_WIN32)
#include <windows.h>
#elif defined(__linux__)
#include <sys/sysinfo.h>
#else
#error "Unsupported platform for getAvailableMemory()"
#endif

size_t getAvailableMemory() {
#if defined(_WIN32)
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    if (!GlobalMemoryStatusEx(&status)) {
        // If unable to retrieve, return 0
        return 0;
    }
    return static_cast<size_t>(status.ullAvailPhys);

#elif defined(__linux__)
    struct sysinfo info;
    if (sysinfo(&info) == 0) {
        return static_cast<size_t>(info.freeram) * info.mem_unit;
    }
    return 0;
#endif
}

std::vector<std::string> readCombinations(const int& dims) {
    std::vector<std::string> comb;

    // Supported dimension range: 2 to 9
    if (dims < 2 || dims > 9) {
        std::cerr << "Error in readCombinations: dimension "
                  << dims << " not supported.\n";
        return comb; // Return empty
    }

    // Recupera il path di *questo* file sorgente, così come è noto a compile-time.
    std::filesystem::__cxx11::path thisSourceFile = std::filesystem::__cxx11::path(__FILE__);
    // Ottiene la cartella che contiene il file sorgente.
    std::filesystem::__cxx11::path thisSourceDir  = thisSourceFile.parent_path();
    // A questo punto, “bin” è fisso rispetto alla cartella del file sorgente.
    std::filesystem::__cxx11::path binDir         = thisSourceDir / "../bin";

    // Costruisci l’array dei nomi di file nella cartella bin
    std::string FileName[] = {
        (binDir / "Comb2D.txt").string(),
        (binDir / "Comb3D.txt").string(),
        (binDir / "Comb4D.txt").string(),
        (binDir / "Comb5D.txt").string(),
        (binDir / "Comb6D.txt").string(),
        (binDir / "Comb7D.txt").string(),
        (binDir / "Comb8D.txt").string(),
        (binDir / "Comb9D.txt").string()
    };
    std::string fileToOpen = FileName[dims - 2];

    FILE* fp = std::fopen(fileToOpen.c_str(), "r");
    if (!fp) {
        std::cerr << "Error: cannot open file " << fileToOpen << std::endl;
        std::exit(1);
    }

    char buf[512];
    // First line should match dims
    if (!std::fgets(buf, 512, fp)) {
        std::cerr << "Error reading from " << fileToOpen << std::endl;
        std::fclose(fp);
        std::exit(1);
    }
    if (std::atoi(buf) != dims) {
        std::cerr << "Error! Dimensions do not match for " << fileToOpen << std::endl;
        std::fclose(fp);
        std::exit(1);
    }

    // Read subsequent lines as binary strings
    while (std::fgets(buf, 512, fp)) {
        char* token = std::strtok(buf, " \n\t");
        if (token) {
            comb.emplace_back(token);
        }
    }
    std::fclose(fp);
    return comb;
}

bool MbrIsValid(const std::vector<std::array<float, 2>>& mbr,
                const std::vector<std::string>& Comb,
                int dims,
                const float queryPlane[])
{
    int numAbove = 0;
    int numBelow = 0;
    long numOfVertices = static_cast<long>(Comb.size());

    // Temporary storage for vertex coords
    std::vector<double> coord(dims, 0.0);
    coord.reserve(dims);

    for (const auto& combination : Comb) {
        float sumVal = 0.0f;
        for (int j = 0; j < dims; ++j) {
            coord[j] = (combination[j] == '0') ? mbr[j][0] : mbr[j][1];
            sumVal += static_cast<float>(coord[j]);
        }
        // Compare sumVal to queryPlane[dims]
        if (sumVal > queryPlane[dims]) {
            ++numAbove;
            // If some vertices are above and some are below => MBR intersects
            if (numAbove > 0 && numBelow > 0) {
                break;
            }
        } else if (sumVal < queryPlane[dims]) {
            ++numBelow;
            if (numAbove > 0 && numBelow > 0) {
                break;
            }
        }
    }

    // If ALL vertices are above => MBR is invalid
    if (numAbove == numOfVertices) {
        return false;
    }
    // If ALL vertices are below => MBR is valid
    if (numBelow == numOfVertices) {
        return true;
    }
    // Otherwise, it intersects => not valid
    return false;
}
