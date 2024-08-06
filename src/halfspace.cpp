//
// Created by leona on 06/08/2024.
//

#include "geom.h"
#include "halfspace.h"
#include <iostream>
#include <sstream>
#include <cassert>
#include <limits>

Point::Point(const std::vector<double>& coord, int id) : id(id), coord(coord), dims(coord.size()) {}

HalfLine::HalfLine(const Point& pnt) : pnt(pnt), dims(2), arr(Arrangement::AUGMENTED) {
    m = pnt.coord[0] - pnt.coord[1];
    q = pnt.coord[1];
}

double HalfLine::get_y(double x) const {
    return m * x + q;
}

HalfSpace::HalfSpace(const Point& pnt, const Eigen::VectorXd& coeff, double known)
    : pnt(pnt), coeff(coeff), known(known), dims(coeff.size()), arr(Arrangement::AUGMENTED) {}


Point find_halflines_intersection(const HalfLine& r, const HalfLine& s) {
    if (r.m == s.m) {
        return Point({std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()});
    } else {
        double x = (s.q - r.q) / (r.m - s.m);
        return Point({x, r.get_y(x)});
    }
}

Position find_pointhalfspace_position(const Point& point, const HalfSpace& halfspace) {
    double val = halfspace.coeff.dot(Eigen::Map<const Eigen::VectorXd>(point.coord.data(), point.dims));

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
    Eigen::VectorXd p_i = Eigen::Map<const Eigen::VectorXd>(p.coord.data(), p.coord.size() - 1);

    for (const auto& r : records) {
        double r_d = r.coord.back();
        Eigen::VectorXd r_i = Eigen::Map<const Eigen::VectorXd>(r.coord.data(), r.coord.size() - 1);

        // less-than form
        // s(r) <= s(p)
        Eigen::VectorXd coeff = r_i - Eigen::VectorXd::Constant(r_i.size(), r_d) - p_i + Eigen::VectorXd::Constant(p_i.size(), p_d);
        halfspaces.emplace_back(r, coeff, p_d - r_d);
    }

    return halfspaces;
}