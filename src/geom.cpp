//
// Created by leona on 01/08/2024.
//

#include "geom.h"

std::array<std::vector<std::vector<double>>, 2> genMasks(int dims) {
    std::array<std::vector<std::vector<double>>, 2> masks;

    //TODO

    // Resize the vectors to hold the appropriate number of elements
    masks[0].resize(9, std::vector<double>(dims)); // Point masks
    masks[1].resize(9, std::vector<double>(dims)); // Node masks

    // Example masks initialization
    for (int i = 0; i < 9; ++i) {
        for (int j = 0; j < dims; ++j) {
            masks[0][i][j] = (i == j) ? 1.0 : 0.0; // Identity-like mask for points
            masks[1][i][j] = (i % 2 == 0) ? 1.0 : 0.0; // Alternating mask for nodes
        }
    }

    return masks;
}