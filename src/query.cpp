#include "query.h"
#include <algorithm>  // Per std::remove_if
#include <memory>

void getdominators(const std::vector<std::shared_ptr<Point>>& data, const Point& p, std::vector<std::shared_ptr<Point>>& dominators) {
    dominators.clear();  // Pulisci il vector

    for (const auto& r : data) {
        bool less_equal = true;
        bool less = false;

        for (int i = 0; i < p.dims; ++i) {
            if (r->coord[i] > p.coord[i]) {
                less_equal = false;
                break;
            }
            if (r->coord[i] < p.coord[i]) {
                less = true;
            }
        }

        if (less_equal && less) {
            dominators.push_back(r);  // Aggiungi dominatori al vector
        }
    }
}

void getdominees(const std::vector<std::shared_ptr<Point>>& data, const Point& p, std::vector<std::shared_ptr<Point>>& dominees) {
    dominees.clear();  // Pulisci il vector

    for (const auto& r : data) {
        bool greater_equal = true;
        bool greater = false;

        for (int i = 0; i < p.dims; ++i) {
            if (r->coord[i] < p.coord[i]) {
                greater_equal = false;
                break;
            }
            if (r->coord[i] > p.coord[i]) {
                greater = true;
            }
        }

        if (greater_equal && greater) {
            dominees.push_back(r);  // Aggiungi dominees al vector
        }
    }
}

void getincomparables(const std::vector<std::shared_ptr<Point>>& data, const Point& p, std::vector<std::shared_ptr<Point>>& incomparables) {
    incomparables.clear();  // Pulisci il vector

    for (const auto& r : data) {
        bool less = false;
        bool greater = false;

        for (int i = 0; i < p.dims; ++i) {
            if (r->coord[i] < p.coord[i]) {
                less = true;
            }
            if (r->coord[i] > p.coord[i]) {
                greater = true;
            }
        }

        if (less && greater) {
            incomparables.push_back(r);  // Aggiungi incomparables al vector
        }
    }
}

void getskyline(const std::vector<std::shared_ptr<Point>>& data, std::vector<std::shared_ptr<Point>>& skyline) {
    skyline.clear();  // Pulisci il vector

    auto dominates = [](const std::shared_ptr<Point>& p, const std::shared_ptr<Point>& r) {
        bool less_equal = true;
        bool less = false;

        for (int i = 0; i < p->dims; ++i) {
            if (p->coord[i] > r->coord[i]) {
                less_equal = false;
                break;
            }
            if (p->coord[i] < r->coord[i]) {
                less = true;
            }
        }

        return less_equal && less;
    };

    for (const auto& pnt : data) {
        bool dominated = false;

        for (const auto& s : skyline) {
            if (dominates(s, pnt)) {
                dominated = true;
                break;
            }
        }

        if (!dominated) {
            // Rimuovi i punti nel skyline che sono dominati da pnt
            skyline.erase(std::remove_if(skyline.begin(), skyline.end(), [&](const std::shared_ptr<Point>& s) {
                return dominates(pnt, s);
            }), skyline.end());

            // Aggiungi pnt allo skyline
            skyline.push_back(pnt);
        }
    }
}
