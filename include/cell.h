//
// Created by leona on 06/08/2024.
//

#ifndef CELL_H
#define CELL_H

#include <algorithm>
#include <string>
#include <utility> // For std::pair
#include <vector>
#include "geom.h"
#include "halfspace.h"
#include "qnode.h"

class Cell {
public:
    Cell(int order, const std::string& mask, const std::vector<long>& covered, const std::vector<long>& halfspaces, const std::vector<std::array<double, 2>>& leaf_mbr, const Point& feasible_pnt);

    bool issingular() const;

    int order;
    std::string mask;
    std::vector<long> covered;
    std::vector<long> halfspaces;
    std::vector<std::array<double, 2>> leaf_mbr;
    Point feasible_pnt;
};

struct Interval {
    // La retta di “supporto” che ha generato questo intervallo
    std::shared_ptr<HalfLine> halfline;

    // L’intervallo [range.first, range.second] sull’asse x
    // dove l'ordine rimane costante
    std::pair<double, double> range;

    // Se vero, allora al termine di questo intervallo
    // la halfline “esce” dalla copertura. Altrimenti “entra”.
    bool coversleft;

    // Quante halflines coprono quell’intervallo
    int order;

    // Tutte le halflines che coprono l’intervallo
    std::vector<std::shared_ptr<HalfLine>> covered;

    Interval(std::shared_ptr<HalfLine> hl,
             std::pair<double, double> r,
             bool c)
        : halfline(std::move(hl)), range(r), coversleft(c), order(0) {}

    bool issingular() const {
        return std::all_of(covered.begin(), covered.end(),
    [](const std::shared_ptr<HalfLine>& hl) { return hl->arr == Arrangement::SINGULAR; });

    }
};

std::vector<std::string> genhammingstrings(int strlen, int weight);
std::vector<Cell> searchmincells_lp(const QNode& leaf, const std::vector<std::string>& hamstrings);

#endif // CELL_H
