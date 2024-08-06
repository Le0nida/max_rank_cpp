//
// Created by leona on 01/08/2024.
//

#ifndef GEOM_H
#define GEOM_H

#include <vector>

class Point {
public:
    Point(const std::vector<double>& coord, int id = -1);
    int id;
    std::vector<double> coord;
    int dims;
};

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> genmasks(int dims);

#endif // GEOM_H
