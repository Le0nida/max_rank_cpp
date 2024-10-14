//
// Created by leona on 01/08/2024.
//

#ifndef HALFSPACE_H
#define HALFSPACE_H

#include "geom.h"
#include <cstdlib> // Per malloc e free
#include <memory>
#include <vector>


enum Position {
    POS_IN = 1,
    POS_OUT = -1,
    POS_ON = 0
};

enum Arrangement {
    SINGULAR = 0,
    AUGMENTED = 1
};

class HalfLine {
public:
    explicit HalfLine(const Point& pnt);
    [[nodiscard]] double get_y(double x) const;

    Point pnt;
    double m;
    double q;
    Arrangement arr;
    int dims;
};


class HalfSpace {
public:
    long int pntID;
    std::vector<double> coeff;
    double known;
    Arrangement arr;
    int dims;

    HalfSpace(long int pntID, const std::vector<double>& coeff, double known, int dims)
        : pntID(pntID), coeff(coeff), known(known), arr(AUGMENTED), dims(dims) {}

    // Equality operator
    bool operator==(const HalfSpace& other) const;
};

// Generate halfspaces from a point and a set of records
std::vector<std::shared_ptr<HalfSpace>> genhalfspaces(const Point& p, const std::vector<std::shared_ptr<Point>>& records);

// Other function declarations
Position find_pointhalfspace_position(const Point& point, const HalfSpace& halfspace);

#endif // HALFSPACE_H
