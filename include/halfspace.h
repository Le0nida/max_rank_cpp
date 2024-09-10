//
// Created by leona on 01/08/2024.
//

#ifndef HALFSPACE_H
#define HALFSPACE_H

#include <vector>
#include "geom.h"

enum Position {
    IN = 1,
    OUT = -1,
    ON = 0
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
    HalfSpace(const Point& pnt, const std::vector<double>& coeff, double known);
    HalfSpace();

    Point pnt;
    std::vector<double> coeff;
    double known;
    Arrangement arr;
    int dims;

    bool operator==(const HalfSpace& other) const;

    // Serializza l'oggetto HalfSpace su disco
    void saveToDisk(std::ofstream& out) const {
        pnt.saveToDisk(out);  // Salva il punto associato

        // Serializza i coefficienti
        size_t coeffSize = coeff.size();
        out.write(reinterpret_cast<const char*>(&coeffSize), sizeof(coeffSize));
        out.write(reinterpret_cast<const char*>(coeff.data()), coeffSize * sizeof(double));

        // Serializza altri membri
        out.write(reinterpret_cast<const char*>(&known), sizeof(known));
        out.write(reinterpret_cast<const char*>(&arr), sizeof(arr));
        out.write(reinterpret_cast<const char*>(&dims), sizeof(dims));
    }

    // Carica l'oggetto HalfSpace da disco
    void loadFromDisk(std::ifstream& in) {
        pnt.loadFromDisk(in);  // Carica il punto associato

        // Carica i coefficienti
        size_t coeffSize;
        in.read(reinterpret_cast<char*>(&coeffSize), sizeof(coeffSize));
        coeff.resize(coeffSize);
        in.read(reinterpret_cast<char*>(coeff.data()), coeffSize * sizeof(double));

        // Carica altri membri
        in.read(reinterpret_cast<char*>(&known), sizeof(known));
        in.read(reinterpret_cast<char*>(&arr), sizeof(arr));
        in.read(reinterpret_cast<char*>(&dims), sizeof(dims));
    }
};

Point find_halflines_intersection(const HalfLine& r, const HalfLine& s);
Position find_pointhalfspace_position(const Point& point, const HalfSpace& halfspace);
std::vector<HalfSpace> genhalfspaces(const Point& p, const std::vector<Point>& records);

#endif //HALFSPACE_H
