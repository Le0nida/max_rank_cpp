//
// Created by leona on 06/08/2024.
//

#include "halfspace.h"
#include <cstring> // Per memcpy
#include <cmath>
#include <map>
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

HalfSpace::HalfSpace(long int pntID, double* coeff, double known, int dims)
    : pntID(pntID), known(known), arr(AUGMENTED), dims(dims)
{
    this->coeff = (double*)malloc(dims * sizeof(double));
    memcpy(this->coeff, coeff, dims * sizeof(double));
}

HalfSpace::~HalfSpace() {
    if (coeff) {
        free(coeff);
        coeff = nullptr;
    }
}

bool HalfSpace::operator==(const HalfSpace& other) const {
    if (pntID != other.pntID || known != other.known || arr != other.arr || dims != other.dims)
        return false;

    for (int i = 0; i < dims; ++i) {
        if (coeff[i] != other.coeff[i])
            return false;
    }
    return true;
}

// Generate halfspaces from a point and a set of records
HalfSpace** genhalfspaces(const Point& p, Point** records, Point** old_records, int numRecords, int numOldRecords, int& numHalfSpaces, std::vector<HalfSpace *>& halfspacesToInsert) {
    int dims = p.dims - 1;
    double p_d = p.coord[p.dims - 1];  // Last coordinate of p
    double* p_i = p.coord;  // Coordinates of p, excluding the last

    // Inizializza la capacità iniziale per l'array halfspaces
    int capacity = numRecords;
    HalfSpace** halfspaces = (HalfSpace**)malloc(capacity * sizeof(HalfSpace*));
    numHalfSpaces = 0; // Counter per i halfspaces effettivamente generati

    for (int idx = 0; idx < numRecords; ++idx) {
        Point* r = records[idx];
        bool found = false;

        if (HalfSpaces.find(r->id) != HalfSpaces.end())
        {
            continue;
        }

        // Calcola i coefficienti dell'halfspace
        double r_d = r->coord[r->dims - 1];
        double* r_i = r->coord;

        double* coeff = (double*)malloc(dims * sizeof(double));
        for (int i = 0; i < dims; ++i) {
            coeff[i] = (r_i[i] - r_d) - (p_i[i] - p_d);
        }

        // Crea un nuovo halfspace con l'id del record corrente
        long int id = r->id;

        auto* hs = new HalfSpace(id, coeff, p_d - r_d, dims);
        halfspacesToInsert.push_back(hs);
        HalfSpaces.emplace(id);
        // Dealloca il coeff poiché è stato copiato nella struttura HalfSpace
        free(coeff);

        // Aggiungi l'halfspace al vettore
        if (numHalfSpaces >= capacity) {
            capacity *= 2;  // Ridimensiona se necessario
            halfspaces = (HalfSpace**)realloc(halfspaces, capacity * sizeof(HalfSpace*));
        }

        halfspaces[numHalfSpaces++] = hs;  // Incrementa il contatore
    }

    // Ridimensiona l'array di halfspaces per rispecchiare il numero esatto di halfspaces creati
    if (numHalfSpaces < capacity) {
        halfspaces = (HalfSpace**)realloc(halfspaces, numHalfSpaces * sizeof(HalfSpace*));
    }

    return halfspaces;
}


Position find_pointhalfspace_position(const Point& point, const HalfSpace& halfspace) {
    double val = 0.0;
    for (int i = 0; i < halfspace.dims; ++i) {
        val += halfspace.coeff[i] * point.coord[i];
    }

    if (val < halfspace.known) {
        return POS_IN;
    } else if (val > halfspace.known) {
        return POS_OUT;
    } else {
        return POS_ON;
    }
}