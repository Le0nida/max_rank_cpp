//
// Created by leona on 01/08/2024.
//

#ifndef GEOM_H
#define GEOM_H

#include <vector>

class Point {
public:
    int id;
    std::vector<double> coord;
    size_t dims;


    explicit Point(const std::vector<double>& coord, const int id = -1)
        : id(id), coord(coord), dims(coord.size()) {}


    // Operatore di uguaglianza per confrontare due punti
    bool operator==(const Point& other) const;
};

// Funzione per generare le maschere
void genmasks(int dims, double***& pts_mask, int& num_pts, double***& nds_mask, int& num_nds);

#endif // GEOM_H

