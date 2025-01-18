#include "geom.h"
#include <array>
#include <bitset>
#include <numeric>

Point::Point(const std::vector<double>& coord, int id)
    : id(id),
      coord(coord),
      dims((int)coord.size())
{
}

std::pair<std::vector<std::vector<double>>,
          std::vector<std::vector<double>>> genmasks(int dims)
{
    // Example function to generate two sets of N-d masks
    std::vector<double> incr(dims, 0.5);
    std::vector<std::vector<double>> pts(1, std::vector<double>(dims, 0.5));

    // Generate points offset by +/- 0.5
    for (int d = 0; d < dims; ++d) {
        std::vector<std::vector<double>> lower = pts;
        std::vector<std::vector<double>> higher = pts;
        for (auto& p : lower)  p[d] -= incr[d];
        for (auto& p : higher) p[d] += incr[d];
        pts.insert(pts.end(), lower.begin(), lower.end());
        pts.insert(pts.end(), higher.begin(), higher.end());
    }

    // Build masks
    std::vector<std::vector<double>> pts_mask(pts.size(), std::vector<double>(dims));
    for (size_t i = 0; i < pts.size(); ++i) {
        for (int d = 0; d < dims; ++d) {
            pts_mask[i][d] = (pts[i][d] - incr[d]) / incr[d];
        }
    }

    // Prepare MBR array
    int total = (1 << dims);
    std::vector<std::vector<std::array<double, 2>>> mbr(total, std::vector<std::array<double, 2>>(dims));
    for (int quad = 0; quad < total; ++quad) {
        std::bitset<32> bits(quad);
        for (int d = 0; d < dims; ++d) {
            double minval = (bits[d] ? 0.5 : 0.0);
            double maxval = (bits[d] ? 1.0 : 0.5);
            mbr[quad][d] = { minval, maxval };
        }
    }

    // Build nds_mask
    std::vector<std::vector<double>> nds_mask(pts.size(),
                                              std::vector<double>(total, 0));
    for (size_t p = 0; p < pts.size(); ++p) {
        for (int n = 0; n < total; ++n) {
            bool match = true;
            for (int d = 0; d < dims; ++d) {
                double val = pts[p][d];
                if (val != mbr[n][d][0] && val != mbr[n][d][1]) {
                    match = false;
                    break;
                }
            }
            if (match) {
                nds_mask[p][n] = 1;
            }
        }
    }

    return { pts_mask, nds_mask };
}
