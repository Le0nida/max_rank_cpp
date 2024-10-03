//
// Created by leona on 06/08/2024.
//

#ifndef CELL_H
#define CELL_H

#include <string>
#include <utility>
#include <vector>
#include "geom.h"
#include "halfspace.h"
#include "qnode.h"

class Context; // Forward declaration

class Cell {
public:
    Cell(int order, const std::string& mask, const std::vector<long>& covered,
         const std::vector<long>& halfspaces, const std::vector<std::array<double, 2>>& leaf_mbr,
         const Point& feasible_pnt);

    bool issingular(Context& ctx) const;

    int order;
    std::string mask;
    std::vector<long> covered;
    std::vector<long> halfspaces;
    std::vector<std::array<double, 2>> leaf_mbr;
    Point feasible_pnt;
};

std::vector<std::string> genhammingstrings(int strlen, int weight);
std::vector<Cell> searchmincells_lp(Context& ctx, const QNode& leaf, const std::vector<std::string>& hamstrings);

#endif // CELL_H