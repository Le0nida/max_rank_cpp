#include "query.h"
#include <algorithm>
#include <cmath>

// ------------------------------------------------------------------------
// Existing helper functions (unchanged)
// ------------------------------------------------------------------------
std::vector<Point> getdominators(const std::vector<Point>& data, const Point& p)
{
    std::vector<Point> dominators;
    dominators.reserve(data.size() / 2);
    for (const auto& r : data) {
        bool less_equal    = true;
        bool strictly_less = false;
        for (int i = 0; i < p.dims; ++i) {
            if (r.coord[i] > p.coord[i]) {
                less_equal = false;
                break;
            }
            if (r.coord[i] < p.coord[i]) {
                strictly_less = true;
            }
        }
        if (less_equal && strictly_less) {
            dominators.push_back(r);
        }
    }
    return dominators;
}

std::vector<Point> getdominees(const std::vector<Point>& data, const Point& p)
{
    std::vector<Point> dominees;
    dominees.reserve(data.size() / 2);
    for (const auto& r : data) {
        bool greater_equal = true;
        bool strictly_greater = false;
        for (int i = 0; i < p.dims; ++i) {
            if (r.coord[i] < p.coord[i]) {
                greater_equal = false;
                break;
            }
            if (r.coord[i] > p.coord[i]) {
                strictly_greater = true;
            }
        }
        if (greater_equal && strictly_greater) {
            dominees.push_back(r);
        }
    }
    return dominees;
}

std::vector<Point> getincomparables(const std::vector<Point>& data, const Point& p)
{
    std::vector<Point> incomp;
    incomp.reserve(data.size() / 2);
    for (const auto& r : data) {
        bool less = false;
        bool greater = false;
        for (int i = 0; i < p.dims; ++i) {
            if (r.coord[i] < p.coord[i]) less = true;
            if (r.coord[i] > p.coord[i]) greater = true;
        }
        if (less && greater) {
            incomp.push_back(r);
        }
    }
    return incomp;
}

/**
 * \brief Controlla se p domina strettamente q (p < q in tutte le dim?).
 *        Ritorna true se p domina q.
 */
static bool dominates(const Point& p, const Point& q)
{
    bool strictlyLess = false;
    for (int i = 0; i < p.dims; i++) {
        if (p.coord[i] > q.coord[i]) {
            // p non domina q
            return false;
        }
        if (p.coord[i] < q.coord[i]) {
            strictlyLess = true;
        }
    }
    return strictlyLess;
}

/**
 * \brief Calcola la skyline di data con l'approccio Sort–Filter–Skyline (SFS).
 *        1) Ordina i punti per la somma delle coordinate, in modo crescente
 *        2) Scansiona i punti e mantiene un insieme di skyline incrementale
 */
std::vector<Point> getskyline(const std::vector<Point>& data)
{
    if (data.empty()) return {};

    // 1) Prepara un vettore (punto, sumCoord) e ordina
    std::vector<std::pair<Point, double>> arr;
    arr.reserve(data.size());
    for (const auto& p : data) {
        double sum = 0.0;
        for (double c : p.coord) {
            sum += c;
        }
        arr.emplace_back(p, sum);
    }

    // Ordiniamo in base a sumCoord (crescente)
    std::sort(arr.begin(), arr.end(),
              [](auto& a, auto& b){
                  return a.second < b.second;
              });

    // 2) Filtro incrementale per costruire la skyline
    std::vector<Point> sky;
    sky.reserve(data.size() / 10);

    for (auto &pairp : arr) {
        Point p = pairp.first;
        bool dominated = false;

        // Se p è dominato da uno qualunque in sky, scartiamo p
        for (auto &sk : sky) {
            if (dominates(sk, p)) {
                dominated = true;
                break;
            }
        }
        if (!dominated) {
            // p non è dominato, lo aggiungiamo e
            // rimuoviamo eventuali punti in sky che p domina
            sky.erase(std::remove_if(sky.begin(), sky.end(),
                                     [&](const Point &x){
                                         return dominates(p, x);
                                     }),
                      sky.end());
            sky.push_back(p);
        }
    }

    return sky;
}