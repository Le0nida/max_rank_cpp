//
// Created by leona on 06/08/2024.
//

#include "halfspace.h"
#include <limits>
#include <numeric>
#include "context.h"

HalfLine::HalfLine(const Point& pnt)
    : pnt(pnt), dims(2), arr(Arrangement::AUGMENTED) {
    m = pnt.coord[0] - pnt.coord[1];
    q = pnt.coord[1];
}

double HalfLine::get_y(double x) const {
    return m * x + q;
}

HalfSpace::HalfSpace(long int pntID, const std::vector<double>& coeff, double known)
    : pntID(pntID), coeff(coeff), known(known), dims(coeff.size()), arr(Arrangement::AUGMENTED) {}

HalfSpace::HalfSpace()
    : pntID(-1), coeff(0), known(0), dims(0), arr(Arrangement::AUGMENTED) {}

bool HalfSpace::operator==(const HalfSpace& other) const {
    return pntID == other.pntID && coeff == other.coeff && known == other.known && arr == other.arr && dims == other.dims;
}

void HalfSpace::saveToDisk(std::ofstream& out) const {
    out.write(reinterpret_cast<const char*>(&pntID), sizeof(pntID));

    size_t coeffSize = coeff.size();
    out.write(reinterpret_cast<const char*>(&coeffSize), sizeof(coeffSize));
    out.write(reinterpret_cast<const char*>(coeff.data()), coeffSize * sizeof(double));

    out.write(reinterpret_cast<const char*>(&known), sizeof(known));
    int arr_int = static_cast<int>(arr);
    out.write(reinterpret_cast<const char*>(&arr_int), sizeof(arr_int));
    out.write(reinterpret_cast<const char*>(&dims), sizeof(dims));
}

void HalfSpace::loadFromDisk(std::ifstream& in) {
    in.read(reinterpret_cast<char*>(&pntID), sizeof(pntID));

    size_t coeffSize;
    in.read(reinterpret_cast<char*>(&coeffSize), sizeof(coeffSize));
    coeff.resize(coeffSize);
    in.read(reinterpret_cast<char*>(coeff.data()), coeffSize * sizeof(double));

    in.read(reinterpret_cast<char*>(&known), sizeof(known));
    int arr_int;
    in.read(reinterpret_cast<char*>(&arr_int), sizeof(arr_int));
    arr = static_cast<Arrangement>(arr_int);
    in.read(reinterpret_cast<char*>(&dims), sizeof(dims));
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
        return Position::POS_IN;
    } else if (val > halfspace.known) {
        return Position::POS_OUT;
    } else {
        return Position::POS_ON;
    }
}

std::vector<long> genhalfspaces(Context& ctx, const Point& p, const std::vector<Point>& records) {
    std::vector<long> halfspaceIDs;
    halfspaceIDs.reserve(records.size());

    double p_d = p.coord.back();
    std::vector<double> p_i(p.coord.begin(), p.coord.end() - 1);

    for (const auto& r : records) {
        if (ctx.pointToHalfSpaceCache.find(r) != ctx.pointToHalfSpaceCache.end()) {
            halfspaceIDs.push_back(ctx.pointToHalfSpaceCache[r]);
            continue;
        }

        double r_d = r.coord.back();
        std::vector<double> r_i(r.coord.begin(), r.coord.end() - 1);

        std::vector<double> coeff(r_i.size());
        for (size_t i = 0; i < r_i.size(); ++i) {
            coeff[i] = (r_i[i] - r_d) - (p_i[i] - p_d);
        }

        long int id = r.id;
        auto halfspace = std::make_shared<HalfSpace>(id, coeff, p_d - r_d);

        ctx.halfspaceCache->insert(id, halfspace);
        halfspaceIDs.push_back(id);
        ctx.pointToHalfSpaceCache[r] = id;
    }

    return halfspaceIDs;
}

HalfSpaceCache::HalfSpaceCache(size_t cacheSize) {
    cache.reserve(cacheSize);
}

void HalfSpaceCache::insert(long id, std::shared_ptr<HalfSpace> halfspace) {
    cache[id] = std::move(halfspace);
}

std::shared_ptr<HalfSpace> HalfSpaceCache::get(long id) const {
    auto it = cache.find(id);
    if (it != cache.end()) {
        return it->second;
    }
    return nullptr;
}

bool HalfSpaceCache::contains(long id) const {
    return cache.find(id) != cache.end();
}