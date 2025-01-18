#include "query.h"
#include <algorithm>

/**
 * \brief Gathers all points that strictly dominate p in at least one dimension
 *        and do not exceed p in others.
 */
std::vector<Point> getdominators(const std::vector<Point>& data,
                                 const Point& p)
{
    std::vector<Point> dominators;
    dominators.reserve(data.size() / 2);

    for (const auto& r : data) {
        bool less_equal = true;
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

/**
 * \brief Gathers all points that are strictly dominated by p.
 */
std::vector<Point> getdominees(const std::vector<Point>& data,
                               const Point& p)
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

/**
 * \brief Gathers points that are neither dominated by p nor dominating p.
 */
std::vector<Point> getincomparables(const std::vector<Point>& data,
                                    const Point& p)
{
    std::vector<Point> incomp;
    incomp.reserve(data.size() / 2);

    for (const auto& r : data) {
        bool less = false;
        bool greater = false;

        for (int i = 0; i < p.dims; ++i) {
            if (r.coord[i] < p.coord[i]) {
                less = true;
            }
            if (r.coord[i] > p.coord[i]) {
                greater = true;
            }
        }
        // If both 'less' and 'greater' are true, r is incomparable
        if (less && greater) {
            incomp.push_back(r);
        }
    }
    return incomp;
}

/**
 * \brief Returns a minimal "skyline" of points: not dominated by any other point.
 */
std::vector<Point> getskyline(const std::vector<Point>& data)
{
    auto dominates = [](const Point& p, const Point& r) {
        bool less_equal = true;
        bool strictly_less = false;

        for (int i = 0; i < p.dims; ++i) {
            if (p.coord[i] > r.coord[i]) {
                less_equal = false;
                break;
            }
            if (p.coord[i] < r.coord[i]) {
                strictly_less = true;
            }
        }
        return (less_equal && strictly_less);
    };

    std::vector<Point> window;
    window.reserve(data.size() / 2);

    for (const auto& pnt : data) {
        bool isDominated = false;
        for (const auto& w_pnt : window) {
            if (dominates(w_pnt, pnt)) {
                isDominated = true;
                break;
            }
        }
        if (!isDominated) {
            // Remove window points that are dominated by pnt
            window.erase(std::remove_if(window.begin(), window.end(),
                            [&](const Point& w_pnt) {
                                return dominates(pnt, w_pnt);
                            }),
                         window.end());
            window.push_back(pnt);
        }
    }
    return window;
}
