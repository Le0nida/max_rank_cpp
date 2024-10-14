//
// Created by leona on 01/08/2024.
//

#include "geom.h"
#include <cmath>
#include <cstring> // Per memcpy
#include <bitset>

bool Point::operator==(const Point& other) const {
    if (dims != other.dims) {
        return false;
    }
    for (int i = 0; i < dims; ++i) {
        if (coord[i] != other.coord[i]) {
            return false;
        }
    }
    return true;
}

// Funzione per generare le maschere
void genmasks(int dims, double***& pts_mask, int& num_pts, double***& nds_mask, int& num_nds) {
    int capacity = 1024;
    num_pts = 1;
    pts_mask = (double***)malloc(capacity * sizeof(double**));
    pts_mask[0] = (double**)malloc(dims * sizeof(double*));
    for (int i = 0; i < dims; ++i) {
        pts_mask[0][i] = (double*)malloc(sizeof(double));
        pts_mask[0][i][0] = 0.5;
    }

    double* incr = (double*)malloc(dims * sizeof(double));
    for (int i = 0; i < dims; ++i) {
        incr[i] = 0.5;
    }

    // Genera i punti
    for (int d = 0; d < dims; ++d) {
        int current_num_pts = num_pts;
        for (int p = 0; p < current_num_pts; ++p) {
            // Punto lower
            auto** lower = (double**)malloc(dims * sizeof(double*));
            // Punto higher
            auto** higher = (double**)malloc(dims * sizeof(double*));

            for (int i = 0; i < dims; ++i) {
                lower[i] = (double*)malloc(sizeof(double));
                higher[i] = (double*)malloc(sizeof(double));
                if (i == d) {
                    lower[i][0] = pts_mask[p][i][0] - incr[i];
                    higher[i][0] = pts_mask[p][i][0] + incr[i];
                } else {
                    lower[i][0] = pts_mask[p][i][0];
                    higher[i][0] = pts_mask[p][i][0];
                }
            }

            // Aggiungi lower e higher a pts_mask
            if (num_pts + 2 > capacity) {
                capacity *= 2;
                pts_mask = (double***)realloc(pts_mask, capacity * sizeof(double**));
            }
            pts_mask[num_pts++] = lower;
            pts_mask[num_pts++] = higher;
        }
    }

    // Calcola pts_mask
    for (int p = 0; p < num_pts; ++p) {
        for (int d = 0; d < dims; ++d) {
            pts_mask[p][d][0] = (pts_mask[p][d][0] - incr[d]) / incr[d];
        }
    }

    // Genera mbr
    num_nds = 1 << dims;
    auto**** mbr = (double****)malloc(num_nds * sizeof(double***));
    for (int quad = 0; quad < num_nds; ++quad) {
        mbr[quad] = (double***)malloc(dims * sizeof(double**));
        for (int d = 0; d < dims; ++d) {
            mbr[quad][d] = (double**)malloc(2 * sizeof(double*));
            mbr[quad][d][0] = (double*)malloc(sizeof(double));
            mbr[quad][d][1] = (double*)malloc(sizeof(double));
        }
        for (int d = 0; d < dims; ++d) {
            if (quad & (1 << d)) {
                mbr[quad][d][0][0] = 0.5;
                mbr[quad][d][1][0] = 1.0;
            } else {
                mbr[quad][d][0][0] = 0.0;
                mbr[quad][d][1][0] = 0.5;
            }
        }
    }

    // Calcola nds_mask
    nds_mask = (double***)malloc(num_pts * sizeof(double**));
    for (int p = 0; p < num_pts; ++p) {
        nds_mask[p] = (double**)malloc(num_nds * sizeof(double*));
        for (int n = 0; n < num_nds; ++n) {
            nds_mask[p][n] = (double*)malloc(sizeof(double));
            bool match = true;
            for (int d = 0; d < dims; ++d) {
                double val = pts_mask[p][d][0];
                if (val != mbr[n][d][0][0] && val != mbr[n][d][1][0]) {
                    match = false;
                    break;
                }
            }
            nds_mask[p][n][0] = match ? 1.0 : 0.0;
        }
    }

    // Dealloca mbr
    for (int quad = 0; quad < num_nds; ++quad) {
        for (int d = 0; d < dims; ++d) {
            free(mbr[quad][d][0]);
            free(mbr[quad][d][1]);
            free(mbr[quad][d]);
        }
        free(mbr[quad]);
    }
    free(mbr);
    free(incr);
}
