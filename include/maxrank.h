//
// Created by leona on 06/08/2024.
//

#ifndef MAXRANK_H
#define MAXRANK_H

#include "geom.h"
#include "cell.h"

// Definisci la funzione aa_hd
std::pair<int, std::vector<std::shared_ptr<Cell>>> aa_hd(const std::vector<std::shared_ptr<Point>>& data, const Point& p);

#endif // MAXRANK_H
