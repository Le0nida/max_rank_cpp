//
// Created by leona on 06/08/2024.
//

#ifndef MAXRANK_H
#define MAXRANK_H
#include <string>
#include <vector>

#include "halfspace.h"
#include "qnode.h"

class HalfSpace;

class Cell {
public:
    Cell(int order, const std::string& mask, const std::vector<HalfSpace>& covered, const std::vector<HalfSpace>& halfspaces, const std::vector<std::array<double, 2>>& leaf_mbr, const Point& feasible_pnt);

    bool issingular() const;

    int order;
    std::string mask;
    std::vector<HalfSpace> covered;
    std::vector<HalfSpace> halfspaces;
    std::vector<std::array<double, 2>> leaf_mbr;
    Point feasible_pnt;
};

class Interval {
public:
    Interval(const HalfLine& halfline, const std::pair<double, double>& range, int coversleft);

    bool issingular() const;

    HalfLine halfline;
    std::pair<double, double> range;
    int coversleft;
    std::vector<HalfSpace> covered;
};

std::vector<std::string> genhammingstrings(int strlen, int weight);
std::vector<Cell> searchmincells_lp(const QNode& leaf, const std::vector<std::string>& hamstrings);

#endif //MAXRANK_H
