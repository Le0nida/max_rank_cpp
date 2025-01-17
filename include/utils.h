//
// Created by leona on 13/09/2024.
//

#ifndef UTILS_H
#define UTILS_H


#include <cstddef>
#include <string>
#include <vector>
#include <array>
#include <iostream>
#include <cstring>
#include <iostream>

size_t getAvailableMemory();

/**
 * Legge da file le combinazioni per le dimensioni date (da 2D a 9D).
 * Il file da cui leggere dipende dal parametro dims.
 * Restituisce un vettore di stringhe contenenti le combinazioni binarie.
 */
std::vector<std::string> readCombinations(const int& dims);

/**
 * Verifica se un MBR (insieme di intervalli [min,max] per dimensione)
 * risulti "valido" rispetto a una iperpiano definito da queryPlane.
 *
 * Comb: insieme di combinazioni binarie ("0"/"1") usate per testare i vertici del MBR.
 * dims: numero di dimensioni effettive (senza l'ultima coordinata "fittizia").
 * queryPlane: array di float che descrive il piano di query (es. q_1 + q_2 + ... + q_d < 1).
 */
bool MbrIsValid(const std::vector<std::array<double, 2>>& mbr,
                const std::vector<std::string>& Comb,
                int dims,
                const float queryPlane[]);


#endif //UTILS_H
