#include "query.h"
#include <vector>
#include <algorithm>

std::vector<Point> getdominators(const std::vector<Point>& data, const Point& p) {
    std::vector<Point> dominators;

    for (const auto& r : data) {
        bool less_equal = true;
        bool less = false;

        for (size_t i = 0; i < p.dims; ++i) {
            if (r.coord[i] > p.coord[i]) {
                less_equal = false;
                break;
            }
            if (r.coord[i] < p.coord[i]) {
                less = true;
            }
        }

        if (less_equal && less) {
            dominators.push_back(r);
        }
    }

    return dominators;
}

std::vector<Point> getdominees(const std::vector<Point>& data, const Point& p) {
    std::vector<Point> dominees;

    for (const auto& r : data) {
        bool greater_equal = true;
        bool greater = false;

        for (size_t i = 0; i < p.dims; ++i) {
            if (r.coord[i] < p.coord[i]) {
                greater_equal = false;
                break;
            }
            if (r.coord[i] > p.coord[i]) {
                greater = true;
            }
        }

        if (greater_equal && greater) {
            dominees.push_back(r);
        }
    }

    return dominees;
}

std::vector<Point> getincomparables(const std::vector<Point>& data, const Point& p) {
    std::vector<Point> incomp;

    for (const auto& r : data) {
        bool less = false;
        bool greater = false;

        for (size_t i = 0; i < p.dims; ++i) {
            if (r.coord[i] < p.coord[i]) {
                less = true;
            }
            if (r.coord[i] > p.coord[i]) {
                greater = true;
            }
        }

        if (less && greater) {
            incomp.push_back(r);
        }
    }

    return incomp;
}

std::vector<Point> getskyline(const std::vector<Point>& data) {
    auto dominates = [](const Point& p, const Point& r) {
        bool less_equal = true;
        bool less = false;

        for (size_t i = 0; i < p.dims; ++i) {
            if (p.coord[i] > r.coord[i]) {
                less_equal = false;
                break;
            }
            if (p.coord[i] < r.coord[i]) {
                less = true;
            }
        }

        return less_equal && less;
    };

    std::vector<Point> window;

    for (const auto& pnt : data) {
        bool dominated = false;
        for (const auto& w_pnt : window) {
            if (dominates(w_pnt, pnt)) {
                dominated = true;
                break;
            }
        }

        if (!dominated) {
            window.erase(std::remove_if(window.begin(), window.end(), [&](const Point& w_pnt) { return dominates(pnt, w_pnt); }), window.end());
            window.push_back(pnt);
        }
    }

    return window;
}
