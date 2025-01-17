//
// Created by leona on 06/08/2024.
//

#ifndef MAXRANK_H
#define MAXRANK_H

#include <vector>
#include "geom.h"
#include "cell.h"
#include "query.h"
#include "halfspace.h"
#include "cell.h"
#include "qtree.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <utils.h>

// Default parameters (extern)
extern int limitHamWeight;
extern int maxLevelQTree;
extern int maxCapacityQNode;
extern int maxNoBinStringToCheck;

std::pair<int, std::vector<Cell>> aa_hd(const std::vector<Point>& data, const Point& p);

std::pair<int, std::vector<Interval>> aa_2d(const std::vector<Point>& data, const Point& p);

#endif // MAXRANK_H
