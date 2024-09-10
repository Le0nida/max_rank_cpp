//
// Created by leona on 01/08/2024.
//

#ifndef GEOM_H
#define GEOM_H

#include <fstream>
#include <vector>

class Point {
public:
    Point(const std::vector<double>& coord, int id = -1);
    int id;
    std::vector<double> coord;
    int dims;

    // Serializza l'oggetto Point su disco
    void saveToDisk(std::ofstream& out) const {
        out.write(reinterpret_cast<const char*>(&id), sizeof(id));
        out.write(reinterpret_cast<const char*>(&dims), sizeof(dims));

        // Serializza le coordinate
        size_t coordSize = coord.size();
        out.write(reinterpret_cast<const char*>(&coordSize), sizeof(coordSize));
        out.write(reinterpret_cast<const char*>(coord.data()), coordSize * sizeof(double));
    }

    // Carica l'oggetto Point da disco
    void loadFromDisk(std::ifstream& in) {
        in.read(reinterpret_cast<char*>(&id), sizeof(id));
        in.read(reinterpret_cast<char*>(&dims), sizeof(dims));

        // Carica le coordinate
        size_t coordSize;
        in.read(reinterpret_cast<char*>(&coordSize), sizeof(coordSize));
        coord.resize(coordSize);
        in.read(reinterpret_cast<char*>(coord.data()), coordSize * sizeof(double));
    }
};

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> genmasks(int dims);

#endif // GEOM_H
