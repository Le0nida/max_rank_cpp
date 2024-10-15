//
// Created by leona on 06/08/2024.
//

#ifndef CELL_H
#define CELL_H

#include <map>
#include <set>
#include <utility>
#include <vector>
#include <string>
#include "geom.h"
#include "halfspace.h"
#include "qnode.h"

// Struttura per il risultato del problema di programmazione lineare
struct LinprogResult {
    double* x;
    double fun;
    int status;
    char message[128];
};

class Interval {
public:
    Interval(const HalfLine& halfline, const std::pair<double, double>& range, int coversleft);

    [[nodiscard]] bool issingular() const;

    HalfLine halfline;
    std::pair<double, double> range;
    int coversleft;
    std::vector<HalfSpace> covered;
};

// Classe Cell ottimizzata
class Cell {
public:
    Cell(int order, const std::string& mask, const std::vector<std::shared_ptr<HalfSpace>>& covered,
        const std::vector<std::shared_ptr<HalfSpace>>& halfspaces, const std::vector<std::pair<double, double>>& leaf_mbr,
        int dims, const Point& feasible_pnt);

    ~Cell() = default; // Lasciamo che i vector gestiscano la memoria

    [[nodiscard]] bool issingular() const;

    bool testHalfspacePair(std::shared_ptr<HalfSpace> HS1, long int IdxHS1, std::shared_ptr<HalfSpace> HS2, long int IdxHS2, const std::vector<std::pair<double, double>>& subDataSpace,
                                 std::multimap<int, std::string> &InValidHammingStr);

    bool GenHammingHalfSpaces(char *OutFileName, const int Dimen, std::vector<char> &HammingString,
                                    std::vector<std::shared_ptr<HalfSpace>> HalfSpaceIDs, const std::vector<std::pair<double, double>>& subDataSpace);

    long int optimizedInNodeIntersection(std::vector<std::pair<long, QNode *> > &Leaves,
                                               std::vector<std::set<std::shared_ptr<HalfSpace>> > &minCellHalfSpaces,
                                               std::vector<std::vector<char> > &binaryString);

    int order;
    std::string mask;
    std::vector<std::shared_ptr<HalfSpace>> covered;
    std::vector<std::shared_ptr<HalfSpace>> halfspaces;
    std::vector<std::pair<double, double>> leaf_mbr;
    int dims;
    Point feasible_pnt;
};

// Funzioni per la programmazione lineare
LinprogResult* linprog_highs(const double* c, const double* A_ub, const double* b_ub,
                             const double* bounds, int num_vars, int num_constraints);

void free_linprog_result(LinprogResult* result);

// Funzioni per generare stringhe di Hamming e cercare minimi cell
char** genhammingstrings(int strlen, int weight, int& numStrings);
std::vector<std::shared_ptr<Cell>> searchmincells_lp(const QNode& leaf, char** hamstrings, int numHamstrings);

#endif // CELL_H
