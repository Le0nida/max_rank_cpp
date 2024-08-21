//
// Created by leona on 01/08/2024.
//

#ifndef HALFSPACE_H
#define HALFSPACE_H

#include <vector>
#include "geom.h"

enum Position {
    IN = 1,
    OUT = -1,
    ON = 0
};

enum class Arrangement {
    SINGULAR = 0,
    AUGMENTED = 1
};


class HalfLine {
public:
    HalfLine(const Point& pnt);
    double get_y(double x) const;

    Point pnt;
    double m;
    double q;
    Arrangement arr;
    int dims;
};

class HalfSpace {
public:
    HalfSpace(const Point& pnt, const std::vector<double>& coeff, double known);

    Point pnt;
    std::vector<double> coeff;
    double known;
    Arrangement arr;
    int dims;

    bool operator==(const HalfSpace& other) const;
};

Point find_halflines_intersection(const HalfLine& r, const HalfLine& s);
Position find_pointhalfspace_position(const Point& point, const HalfSpace& halfspace);
std::vector<HalfSpace> genhalfspaces(const Point& p, const std::vector<Point>& records);

#endif //HALFSPACE_H
