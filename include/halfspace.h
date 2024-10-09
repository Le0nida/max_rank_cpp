//
// Created by leona on 01/08/2024.
//

#ifndef HALFSPACE_H
#define HALFSPACE_H

#include "geom.h"
#include <cstdlib> // Per malloc e free

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
    double* coeff; // Coefficients (C-style array)
    double known;
    Arrangement arr;
    int dims;

    // Constructor
    HalfSpace(long int pntID, double* coeff, double known, int dims);

    // Destructor
    ~HalfSpace();

    // Disable copy constructor and assignment
    HalfSpace(const HalfSpace&) = delete;
    HalfSpace& operator=(const HalfSpace&) = delete;

    // Equality operator
    bool operator==(const HalfSpace& other) const;
};

// Generate halfspaces from a point and a set of records
HalfSpace** genhalfspaces(const Point& p, Point** records, Point** old_records, int numRecords, int numOldRecords, int& numHalfSpaces);

// Other function declarations
Position find_pointhalfspace_position(const Point& point, const HalfSpace& halfspace);
Point find_halflines_intersection(const HalfLine& r, const HalfLine& s);

#endif // HALFSPACE_H
