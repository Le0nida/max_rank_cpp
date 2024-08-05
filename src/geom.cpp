//
// Created by leona on 01/08/2024.
//

#include "geom.h"
#include <vector>
#include <array>
#include <numeric>
#include <bitset>
#include <iostream>

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> genmasks(int dims) {
    std::vector<double> incr(dims, 0.5);
    std::vector<std::vector<double>> pts(1, std::vector<double>(dims, 0.5));

    // Generate points
    for (int d = 0; d < dims; ++d) {
        std::vector<std::vector<double>> lower = pts;
        std::vector<std::vector<double>> higher = pts;
        for (auto& p : lower) p[d] -= incr[d];
        for (auto& p : higher) p[d] += incr[d];
        pts.insert(pts.end(), lower.begin(), lower.end());
        pts.insert(pts.end(), higher.begin(), higher.end());
    }

    // Calculate pts_mask
    std::vector<std::vector<double>> pts_mask(pts.size(), std::vector<double>(dims));
    for (size_t i = 0; i < pts.size(); ++i) {
        for (int d = 0; d < dims; ++d) {
            pts_mask[i][d] = (pts[i][d] - incr[d]) / incr[d];
        }
    }

    // Generate mbr
    std::vector<std::vector<std::array<double, 2>>> mbr(1 << dims, std::vector<std::array<double, 2>>(dims));
    for (int quad = 0; quad < (1 << dims); ++quad) {
        std::bitset<32> qbin(quad);
        std::vector<double> child_mindim(dims);
        std::vector<double> child_maxdim(dims);

        for (int d = 0; d < dims; ++d) {
            child_mindim[d] = qbin[d] ? 0.5 : 0.0;
            child_maxdim[d] = qbin[d] ? 1.0 : 0.5;
        }

        for (int d = 0; d < dims; ++d) {
            mbr[quad][d] = {child_mindim[d], child_maxdim[d]};
        }
    }

    // Calculate nds_mask
    std::vector<std::vector<double>> nds_mask(pts.size(), std::vector<double>(1 << dims, 0));
    for (size_t p = 0; p < pts.size(); ++p) {
        for (int n = 0; n < (1 << dims); ++n) {
            bool match = true;
            for (int d = 0; d < dims; ++d) {
                if (pts[p][d] != mbr[n][d][0] && pts[p][d] != mbr[n][d][1]) {
                    match = false;
                    break;
                }
            }
            if (match) {
                nds_mask[p][n] = 1;
            }
        }
    }

    return {pts_mask, nds_mask};
}
