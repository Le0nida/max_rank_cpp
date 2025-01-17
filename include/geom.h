//
// Created by leona on 01/08/2024.
//

#ifndef GEOM_H
#define GEOM_H

#include <fstream>
#include <vector>

class Point {
public:
    Point(const std::vector<double>& coord, int id = -1);
    int id;
    std::vector<double> coord;
    int dims;

    // Operatore di uguaglianza per confrontare due punti
    bool operator==(const Point& other) const {
        return coord == other.coord;  // Confronta i vettori di coordinate
    }
};

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> genmasks(int dims);

#endif // GEOM_H
