//
// Created by leona on 06/08/2024.
//

#include "query.h"
#include <Eigen/Dense>

std::vector<Point> getdominators(const std::vector<Point>& data, const Point& p) {
    std::vector<Point> dominators;

    for (const auto& r : data) {
        Eigen::Array<bool, Eigen::Dynamic, 1> less_equal = Eigen::Map<const Eigen::VectorXd>(r.coord.data(), r.dims).array() <= Eigen::Map<const Eigen::VectorXd>(p.coord.data(), p.dims).array();
        Eigen::Array<bool, Eigen::Dynamic, 1> less = Eigen::Map<const Eigen::VectorXd>(r.coord.data(), r.dims).array() < Eigen::Map<const Eigen::VectorXd>(p.coord.data(), p.dims).array();

        if (less_equal.all() && less.any()) {
            dominators.push_back(r);
        }
    }

    return dominators;
}

std::vector<Point> getdominees(const std::vector<Point>& data, const Point& p) {
    std::vector<Point> dominees;

    for (const auto& r : data) {
        Eigen::Array<bool, Eigen::Dynamic, 1> greater_equal = Eigen::Map<const Eigen::VectorXd>(r.coord.data(), r.dims).array() >= Eigen::Map<const Eigen::VectorXd>(p.coord.data(), p.dims).array();
        Eigen::Array<bool, Eigen::Dynamic, 1> greater = Eigen::Map<const Eigen::VectorXd>(r.coord.data(), r.dims).array() > Eigen::Map<const Eigen::VectorXd>(p.coord.data(), p.dims).array();

        if (greater_equal.all() && greater.any()) {
            dominees.push_back(r);
        }
    }

    return dominees;
}

std::vector<Point> getincomparables(const std::vector<Point>& data, const Point& p) {
    std::vector<Point> incomp;

    for (const auto& r : data) {
        Eigen::Array<bool, Eigen::Dynamic, 1> less = Eigen::Map<const Eigen::VectorXd>(r.coord.data(), r.dims).array() < Eigen::Map<const Eigen::VectorXd>(p.coord.data(), p.dims).array();
        Eigen::Array<bool, Eigen::Dynamic, 1> greater = Eigen::Map<const Eigen::VectorXd>(r.coord.data(), r.dims).array() > Eigen::Map<const Eigen::VectorXd>(p.coord.data(), p.dims).array();

        if (less.any() && greater.any()) {
            incomp.push_back(r);
        }
    }

    return incomp;
}

std::vector<Point> getskyline(const std::vector<Point>& data) {
    auto dominates = [](const Point& p, const Point& r) {
        Eigen::Array<bool, Eigen::Dynamic, 1> less_equal = Eigen::Map<const Eigen::VectorXd>(p.coord.data(), p.dims).array() <= Eigen::Map<const Eigen::VectorXd>(r.coord.data(), r.dims).array();
        Eigen::Array<bool, Eigen::Dynamic, 1> less = Eigen::Map<const Eigen::VectorXd>(p.coord.data(), p.dims).array() < Eigen::Map<const Eigen::VectorXd>(r.coord.data(), r.dims).array();

        return less_equal.all() && less.any();
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
