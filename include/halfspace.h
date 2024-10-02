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

    // Serializza l'oggetto HalfSpace su disco
    void saveToDisk(std::ofstream& out) const {
        out.write(reinterpret_cast<const char*>(&pntID), sizeof(pntID));

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
        in.read(reinterpret_cast<char*>(&pntID), sizeof(pntID));

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
std::vector<long> genhalfspaces(const Point& p, const std::vector<Point>& records);

// Struttura della cache globale per memorizzare gli half-spaces
class HalfSpaceCache {
public:

    explicit HalfSpaceCache(size_t cacheSize) {
        // Prealloca la dimensione della cache, in base al numero di record (dimensione fissa)
        cache.reserve(cacheSize);
    }

    // Inserisce un nuovo halfspace nella cache
    void insert(long id, std::shared_ptr<HalfSpace> halfspace) {
        cache[id] = std::move(halfspace);
    }

    // Restituisce il puntatore all'halfspace dato l'ID
    std::shared_ptr<HalfSpace> get(long id) const {
        auto it = cache.find(id);
        if (it != cache.end()) {
            return it->second;
        }
        return nullptr;  // Se l'ID non è trovato
    }

    // Controlla se un ID esiste nella cache
    bool contains(long id) const {
        return cache.find(id) != cache.end();
    }

private:
    std::unordered_map<long, std::shared_ptr<HalfSpace>> cache;  // Cache ordinata per ID
};

struct PointHash {
    std::size_t operator()(const Point& p) const {
        std::size_t seed = 0;
        for (double coord : p.coord) {
            seed ^= std::hash<double>{}(coord) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};


// Cache globale degli half-spaces con dimensione fissa
extern HalfSpaceCache* halfspaceCache;

static void initializeCache(const size_t cacheSize) {
    if (halfspaceCache == nullptr) {
        // Inizializza la cache con la dimensione fissa
        halfspaceCache = new HalfSpaceCache(cacheSize);
    }
}

// Cache per i punti già convertiti in half-space, usando unordered_map per efficienza
extern std::unordered_map<Point, long, PointHash> pointToHalfSpaceCache;

#endif //HALFSPACE_H
