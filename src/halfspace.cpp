#include "halfspace.h"
#include <cassert>
#include <limits>
#include <memory>
#include <numeric>

/**
 * \brief Constructs a HalfLine from a 2D point.
 */
HalfLine::HalfLine(const Point& pnt)
    : pnt(pnt),
      m(pnt.coord[0] - pnt.coord[1]),
      q(pnt.coord[1]),
      arr(Arrangement::AUGMENTED),
      dims(2)
{
}

double HalfLine::get_y(double x) const {
    return m * x + q;
}

HalfSpace::HalfSpace(long int pntID,
                     const std::vector<double>& coeff,
                     double known)
    : pntID(pntID),
      coeff(coeff),
      known(known),
      arr(Arrangement::AUGMENTED),
      dims((int)coeff.size())
{
}

HalfSpace::HalfSpace()
    : pntID(-1),
      coeff(),
      known(0.0),
      arr(Arrangement::AUGMENTED),
      dims(0)
{
}

bool HalfSpace::operator==(const HalfSpace& other) const {
    return (pntID == other.pntID &&
            coeff == other.coeff &&
            known == other.known &&
            arr == other.arr &&
            dims == other.dims);
}

Point find_halflines_intersection(const HalfLine& r, const HalfLine& s) {
    // If same slope, no intersection in finite plane
    if (r.m == s.m) {
        return Point({ std::numeric_limits<double>::infinity(),
                       std::numeric_limits<double>::infinity() });
    } else {
        double x = (s.q - r.q) / (r.m - s.m);
        return Point({ x, r.get_y(x) });
    }
}

Position find_pointhalfspace_position(const Point& point,
                                      const HalfSpace& halfspace)
{
    double val = std::inner_product(
                     halfspace.coeff.begin(),
                     halfspace.coeff.end(),
                     point.coord.begin(),
                     0.0
                 );

    if (val < halfspace.known) {
        return Position::POS_IN;
    } else if (val > halfspace.known) {
        return Position::POS_OUT;
    } else {
        return Position::POS_ON;
    }
}

std::vector<long> genhalfspaces(const Point& p,
                                const std::vector<Point>& records)
{
    std::vector<long> halfspaceIDs;
    halfspaceIDs.reserve(records.size());

    double p_d = p.coord.back();  // Last coordinate
    std::vector<double> p_i(p.coord.begin(), p.coord.end() - 1);

    for (const auto& r : records) {
        // Check if already cached
        auto it = pointToHalfSpaceCache.find(r);
        if (it != pointToHalfSpaceCache.end()) {
            halfspaceIDs.push_back(it->second);
            // Skip creation if found in cache
            continue;
        }

        double r_d = r.coord.back();
        std::vector<double> r_i(r.coord.begin(), r.coord.end() - 1);

        // Build halfspace coefficients
        std::vector<double> coeff(r_i.size());
        for (size_t i = 0; i < r_i.size(); ++i) {
            coeff[i] = (r_i[i] - r_d) - (p_i[i] - p_d);
        }

        long id = r.id;
        // Create a new HalfSpace
        auto halfspace = std::make_shared<HalfSpace>(id, coeff, p_d - r_d);

        // Insert into cache
        halfspaceCache->insert(id, halfspace);

        // Save mapping
        halfspaceIDs.push_back(id);
        pointToHalfSpaceCache[r] = id;
    }

    return halfspaceIDs;
}
