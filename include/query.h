//
// Created by leona on 06/08/2024.
//

#ifndef QUERY_H
#define QUERY_H

#include <memory>

#include "geom.h"
#include <vector>

void getdominators(const std::vector<std::shared_ptr<Point>>& data, const Point& p, std::vector<std::shared_ptr<Point>>& dominators);
void getdominees(const std::vector<std::shared_ptr<Point>>& data, const Point& p, std::vector<std::shared_ptr<Point>>& dominees);
void getincomparables(const std::vector<std::shared_ptr<Point>>& data, const Point& p, std::vector<std::shared_ptr<Point>>& incomparables);
void getskyline(const std::vector<std::shared_ptr<Point>>& data, std::vector<std::shared_ptr<Point>>& skyline);

#endif // QUERY_H

