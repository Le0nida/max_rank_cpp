//
// Created by leona on 06/08/2024.
//
#ifndef MAXRANK_H
#define MAXRANK_H

#include "geom.h"
#include "cell.h"
#include <vector>
#include <utility>

class Context; // Forward declaration

std::pair<int, std::vector<Cell>> aa_hd(Context& ctx, const std::vector<Point>& data, const Point& p);

#endif // MAXRANK_H