//
// Created by leona on 06/08/2024.
//

#include "geom.h"
#include "halfspace.h"
#include <iostream>
#include <cassert>
#include <limits>
#include <memory>
#include <numeric> // Per std::inner_product


HalfLine::HalfLine(const Point& pnt) : pnt(pnt), dims(2), arr(Arrangement::AUGMENTED) {
    m = pnt.coord[0] - pnt.coord[1];
    q = pnt.coord[1];
}

double HalfLine::get_y(double x) const {
    return m * x + q;
}

HalfSpace::HalfSpace(const long int pntID, const std::vector<double>& coeff, double known)
    : pntID(pntID), coeff(coeff), known(known), dims(coeff.size()), arr(Arrangement::AUGMENTED) {}

HalfSpace::HalfSpace()
    : pntID(-1), coeff(0), known(0), dims(0), arr(Arrangement::AUGMENTED) {}

bool HalfSpace::operator==(const HalfSpace& other) const {
    return pntID == other.pntID && coeff == other.coeff && known == other.known && arr == other.arr && dims == other.dims;
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

std::vector<long> genhalfspaces(const Point& p, const std::vector<Point>& records) {
    // Preallocazione della lista degli halfspaceID in base al numero di records (dimensione fissa)
    std::vector<long> halfspaceIDs;
    halfspaceIDs.reserve(records.size());

    double p_d = p.coord.back();  // Ultima coordinata di p (costante)
    std::vector<double> p_i(p.coord.begin(), p.coord.end() - 1);  // Coordinate rimanenti di p

    for (const auto& r : records) {
        // Verifica se il punto "r" è già stato convertito in half-space
        if (pointToHalfSpaceCache.find(r) != pointToHalfSpaceCache.end()) {
            // Se esiste già, aggiungi solo l'ID dell'halfspace corrispondente
            halfspaceIDs.push_back(pointToHalfSpaceCache[r]);
            // TODO modifica riportando solo i nuovi hs
            continue;  // Salta la creazione dell'halfspace per questo record
        }

        double r_d = r.coord.back();  // Ultima coordinata di r (costante)
        std::vector<double> r_i(r.coord.begin(), r.coord.end() - 1);  // Coordinate rimanenti di r

        // Calcolo dei coefficienti dell'halfspace
        std::vector<double> coeff(r_i.size());
        for (size_t i = 0; i < r_i.size(); ++i) {
            coeff[i] = (r_i[i] - r_d) - (p_i[i] - p_d);
        }

        long int id = r.id;
        // Crea il nuovo halfspace
        auto halfspace = std::make_shared<HalfSpace>(id, coeff, p_d - r_d);

        // Inserisci l'halfspace nella cache globale
        halfspaceCache->insert(id, halfspace);

        // Aggiungi l'ID dell'halfspace alla lista e memorizza il mapping nel punto cache
        halfspaceIDs.push_back(id);
        pointToHalfSpaceCache[r] = id;  // Mappa il punto all'ID dell'halfspace
    }

    // Restituisci la lista degli ID degli halfspaces
    return halfspaceIDs;
}