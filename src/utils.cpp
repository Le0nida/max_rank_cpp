//
// Created by leona on 13/09/2024.
//

#include "utils.h"


std::vector<std::string> readCombinations(const int& dims) {
    std::vector<std::string> comb;

    FILE* fp;
    char* token;
    char m_separator[] = " \n\t";
    char buf[512];

    // Mappatura file da 2D a 9D (puoi aggiungere altri se necessario)
    std::string FileName[] = {
        "../bin/Comb2D.txt", "../bin/Comb3D.txt", "../bin/Comb4D.txt",
        "../bin/Comb5D.txt", "../bin/Comb6D.txt", "../bin/Comb7D.txt",
        "../bin/Comb8D.txt", "../bin/Comb9D.txt"
    };

    if (dims < 2 || dims > 9) {
        std::cerr << "Error in readCombinations: dimensioni non supportate (" << dims << ").\n";
        return comb;  // Ritorna vuoto
    }

    fp = fopen(FileName[dims - 2].c_str(), "r");
    if (fp == nullptr) {
        std::cout << "error in file open!" << std::endl;
        std::exit(0);
    }

    // Prima riga: controlla che il valore letto corrisponda a dims
    if (fgets(buf, 512, fp) == nullptr) {
        std::cout << "Error reading file " << FileName[dims - 2] << std::endl;
        fclose(fp);
        std::exit(0);
    }
    if (std::atoi(buf) != dims) {
        std::cout << "Error! Dimensions are not equal!" << std::endl;
        fclose(fp);
        std::exit(0);
    }

    // Lettura delle combinazioni
    while (fgets(buf, 512, fp) != nullptr) {
        token = std::strtok(buf, m_separator);
        if (token) {
            comb.emplace_back(token);
        }
    }
    fclose(fp);

    return comb;
}

bool MbrIsValid(const std::vector<std::array<double, 2>>& mbr,
                const std::vector<std::string>& Comb,
                int dims,
                const float queryPlane[])
{
    int numAbove = 0;
    int numBelow = 0;
    const long int numOfVertices = static_cast<long int>(Comb.size());

    // Per evitare continue riallocazioni
    std::vector<double> coord(dims, 0.0);
    coord.reserve(dims);

    for (const auto &combination : Comb) {
        float sum = 0.0f;
        for (int j = 0; j < dims; ++j) {
            coord[j] = (combination[j] == '0') ? mbr[j][0] : mbr[j][1];
            sum += static_cast<float>(coord[j]);
        }

        // Confronta sum con queryPlane[dims]
        if (sum > queryPlane[dims]) {
            ++numAbove;
            if (numAbove > 0 && numBelow > 0) {
                // Early exit: l'MBR è intersecato
                break;
            }
        }
        else if (sum < queryPlane[dims]) {
            ++numBelow;
            if (numAbove > 0 && numBelow > 0) {
                // Early exit: l'MBR è intersecato
                break;
            }
        }
    }

    // Restituisci il risultato in base ai conteggi
    if (numAbove == numOfVertices) return false;   // tutto sopra
    if (numBelow == numOfVertices) return true;    // tutto sotto
    return false; // altrimenti intersecato
}

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