//
// Created by leona on 06/08/2024.
//

#ifndef MAXRANK_H
#define MAXRANK_H

#include "geom.h"
#include "cell.h"

// Definisci la funzione aa_hd
std::pair<int, Cell**> aa_hd(Point** data, int data_size, const Point& p, int& numMinCellsToReturn);

#endif // MAXRANK_H
