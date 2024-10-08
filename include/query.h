//
// Created by leona on 06/08/2024.
//

#ifndef QUERY_H
#define QUERY_H

#include "geom.h"

void getdominators(Point** data, int data_size, const Point& p, Point*** dominators, int& num_dominators);
void getdominees(Point** data, int data_size, const Point& p, Point*** dominees, int& num_dominees);
void getincomparables(Point** data, int data_size, const Point& p, Point*** incomparables, int& num_incomparables);
void getskyline(Point** data, int data_size, Point*** skyline, int& num_skyline);

#endif // QUERY_H
