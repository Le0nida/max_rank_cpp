//
// Created by leona on 01/08/2024.
//

#ifndef GEOM_H
#define GEOM_H

#include <fstream>
#include <vector>
#include <functional>

class Point {
public:
    Point(const std::vector<double>& coord, int id = -1);
    int id;
    std::vector<double> coord;
    int dims;

    bool operator==(const Point& other) const {
        return coord == other.coord;
    }

    // Serialization methods
    void saveToDisk(std::ofstream& out) const;
    void loadFromDisk(std::ifstream& in);
};

struct PointHash {
    std::size_t operator()(const Point& p) const;
};

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> genmasks(int dims);

#endif // GEOM_H
