//
// Created by leona on 01/08/2024.
//

#ifndef HALFSPACE_H
#define HALFSPACE_H

#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>
#include "geom.h"

class Context; // Forward declaration

enum Position {
    POS_IN = 1,
    POS_OUT = -1,
    POS_ON = 0
};

enum class Arrangement {
    SINGULAR = 0,
    AUGMENTED = 1
};

class HalfLine {
public:
    HalfLine(const Point& pnt);
    double get_y(double x) const;

    Point pnt;
    double m;
    double q;
    Arrangement arr;
    int dims;
};

class HalfSpace {
public:
    HalfSpace(long int pntID, const std::vector<double>& coeff, double known);
    HalfSpace();

    long int pntID;
    std::vector<double> coeff;
    double known;
    Arrangement arr;
    int dims;

    bool operator==(const HalfSpace& other) const;

    // Serialization methods
    void saveToDisk(std::ofstream& out) const;
    void loadFromDisk(std::ifstream& in);
};

Point find_halflines_intersection(const HalfLine& r, const HalfLine& s);
Position find_pointhalfspace_position(const Point& point, const HalfSpace& halfspace);
std::vector<long> genhalfspaces(Context& ctx, const Point& p, const std::vector<Point>& records);

// HalfSpaceCache class
class HalfSpaceCache {
public:
    explicit HalfSpaceCache(size_t cacheSize);

    void insert(long id, std::shared_ptr<HalfSpace> halfspace);
    std::shared_ptr<HalfSpace> get(long id) const;
    bool contains(long id) const;

private:
    std::unordered_map<long, std::shared_ptr<HalfSpace>> cache;
};

#endif // HALFSPACE_H
