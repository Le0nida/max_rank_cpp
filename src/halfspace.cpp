//
// Created by leona on 06/08/2024.
//

#include "geom.h"
#include "halfspace.h"
#include <iostream>
#include <sstream>
#include <cassert>
#include <limits>
#include <numeric> // Per std::inner_product


HalfLine::HalfLine(const Point& pnt) : pnt(pnt), dims(2), arr(Arrangement::AUGMENTED) {
    m = pnt.coord[0] - pnt.coord[1];
    q = pnt.coord[1];
}

double HalfLine::get_y(double x) const {
    return m * x + q;
}

HalfSpace::HalfSpace(const Point& pnt, const std::vector<double>& coeff, double known)
    : pnt(pnt), coeff(coeff), known(known), dims(coeff.size()), arr(Arrangement::AUGMENTED) {}

HalfSpace::HalfSpace()
    : pnt({}), coeff(0), known(0), dims(0), arr(Arrangement::AUGMENTED) {}

bool HalfSpace::operator==(const HalfSpace& other) const {
    return pnt.coord == other.pnt.coord && coeff == other.coeff && known == other.known && arr == other.arr && dims == other.dims;
}

Point find_halflines_intersection(const HalfLine& r, const HalfLine& s) {
    if (r.m == s.m) {
        return Point({std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()});
    } else {
        double x = (s.q - r.q) / (r.m - s.m);
        return Point({x, r.get_y(x)});
    }
}

Position find_pointhalfspace_position(const Point& point, const HalfSpace& halfspace) {
    double val = std::inner_product(halfspace.coeff.begin(), halfspace.coeff.end(), point.coord.begin(), 0.0);

    if (val < halfspace.known) {
        return Position::IN;
    } else if (val > halfspace.known) {
        return Position::OUT;
    } else {
        return Position::ON;
    }
}

std::vector<HalfSpace> genhalfspaces(const Point& p, const std::vector<Point>& records) {
    std::vector<HalfSpace> halfspaces;
    double p_d = p.coord.back();
    std::vector<double> p_i(p.coord.begin(), p.coord.end() - 1);

    for (const auto& r : records) {
        double r_d = r.coord.back();
        std::vector<double> r_i(r.coord.begin(), r.coord.end() - 1);

        std::vector<double> coeff(r_i.size());
        for (size_t i = 0; i < r_i.size(); ++i) {
            coeff[i] = (r_i[i] - r_d) - (p_i[i] - p_d);
        }

        halfspaces.emplace_back(r, coeff, p_d - r_d);
    }

    return halfspaces;
}