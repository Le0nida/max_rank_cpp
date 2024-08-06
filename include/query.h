//
// Created by leona on 06/08/2024.
//

#ifndef QUERY_H
#define QUERY_H

#include <vector>
#include "geom.h"

std::vector<Point> getdominators(const std::vector<Point>& data, const Point& p);
std::vector<Point> getdominees(const std::vector<Point>& data, const Point& p);
std::vector<Point> getincomparables(const std::vector<Point>& data, const Point& p);
std::vector<Point> getskyline(const std::vector<Point>& data);

#endif // QUERY_H
