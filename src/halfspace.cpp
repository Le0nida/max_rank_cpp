//
// Created by leona on 06/08/2024.
//

#include "halfspace.h"
#include <cstring> // Per memcpy
#include <cmath>
#include <map>
#include <memory>
#include <numeric> // Per inner_product
#include <unordered_set>
#include <vector>

extern std::unordered_set<long int> HalfSpaces;

HalfLine::HalfLine(const Point& pnt) : pnt(pnt), dims(2), arr(Arrangement::AUGMENTED) {
    m = pnt.coord[0] - pnt.coord[1];
    q = pnt.coord[1];
}

double HalfLine::get_y(double x) const {
    return m * x + q;
}

bool HalfSpace::operator==(const HalfSpace& other) const {
    return pntID == other.pntID && known == other.known && arr == other.arr && coeff == other.coeff;
}

// Generate halfspaces from a point and a set of records
std::vector<std::shared_ptr<HalfSpace>> genhalfspaces(const Point& p, const std::vector<std::shared_ptr<Point>>& records) {
    int dims = p.dims - 1;
    double p_d = p.coord.back();  // Last coordinate of p
    std::vector<std::shared_ptr<HalfSpace>> halfspaces;

    for (const auto& r : records) {
        if (HalfSpaces.find(r->id) != HalfSpaces.end()) continue;

        double r_d = r->coord.back();  // Last coordinate of r
        std::vector<double> coeff(dims);

        for (int i = 0; i < dims; ++i) {
            coeff[i] = (r->coord[i] - r_d) - (p.coord[i] - p_d);
        }

        auto hs = std::make_shared<HalfSpace>(r->id, coeff, p_d - r_d, dims);
        halfspaces.push_back(hs);
        HalfSpaces.insert(r->id);
    }

    return halfspaces;
}

Position find_pointhalfspace_position(const Point& point, const HalfSpace& halfspace) {
    double val = std::inner_product(halfspace.coeff.begin(), halfspace.coeff.end(), point.coord.begin(), 0.0);

    if (val < halfspace.known) return POS_IN;
    if (val > halfspace.known) return POS_OUT;
    return POS_ON;
}