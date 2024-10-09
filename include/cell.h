//
// Created by leona on 06/08/2024.
//

#ifndef CELL_H
#define CELL_H

#include <utility>
#include <vector>

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

    bool issingular() const;

    HalfLine halfline;
    std::pair<double, double> range;
    int coversleft;
    std::vector<HalfSpace> covered;
};

// Classe Cell ottimizzata
class Cell {
public:
    // Costruttore
    Cell(int order, const char* mask, HalfSpace** covered, int numCovered,
         std::vector<long int> halfspaces, double** leaf_mbr, int dims,
         const Point& feasible_pnt);

    // Distruttore
    ~Cell();

    // Funzione per verificare se il cell è singolare
    bool issingular() const;

    int order;
    char* mask;
    HalfSpace** covered;
    int numCovered;
    std::vector<long int> halfspaces;
    double** leaf_mbr;
    int dims;
    Point feasible_pnt;
};

// Funzioni per la programmazione lineare
LinprogResult* linprog_highs(const double* c, const double* A_ub, const double* b_ub,
                             const double* bounds, int num_vars, int num_constraints);

void free_linprog_result(LinprogResult* result);

// Funzioni per generare stringhe di Hamming e cercare minimi cell
char** genhammingstrings(int strlen, int weight, int& numStrings);
Cell** searchmincells_lp(const QNode& leaf, char** hamstrings, int numHamstrings, int& numCells);

#endif // CELL_H
