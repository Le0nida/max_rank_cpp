//
// Created by leona on 01/08/2024.
//

#include "geom.h"

std::vector<std::array<std::vector<double>, 2>> genMasks(int dims) {
    std::vector<std::array<std::vector<double>, 2>> masks(dims);
    // TODO
    for (int i = 0; i < dims; ++i) {
        masks[i][0] = {1.0, 0.0}; // Example mask
        masks[i][1] = {0.0, 1.0}; // Example mask
    }
    return masks;
}