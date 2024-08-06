//
// Created by leona on 06/08/2024.
//

#ifndef MAXRANK_H
#define MAXRANK_H

#include <vector>
#include "geom.h"
#include "qtree.h"
#include "cell.h"

class Cell;
std::pair<int, std::vector<Cell>> aa_hd(const std::vector<Point>& data, const Point& p);

#endif // MAXRANK_H
