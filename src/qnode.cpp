//
// Created by leona on 01/08/2024.
//

#include "qnode.h"
#include <cmath>
#include <map>

// Variabili esterne
int normalizedMax = 1;
int maxCapacity = 10; // Capacità massima del nodo
int globalNodeID = 0;
int maxLevel = 6;
extern int numOfSubdivisions;

QNode::QNode(QNode* parent, const std::vector<std::pair<double, double>>& mbr, int dims)
    : parent(parent), mbr(mbr), norm(true), leaf(true), order(0), dims(dims) {
    nodeID = globalNodeID++;
    level = parent == nullptr ? 0 : parent->level + 1;
}

void QNode::setOrder() {
    size_t localOrder = covered.size();
    QNode* ref = parent;

    while (ref != nullptr) {
        localOrder += ref->covered.size();
        ref = ref->parent;
    }

    order = localOrder;
}

PositionHS QNode::MbrVersusHalfSpace(const std::vector<double>& hs_coeff, const double hs_known) {
    double minVal = 0.0;
    double maxVal = 0.0;

    // Primo ciclo: gestisci coefficienti positivi
    for (int i = 0; i < dims; ++i) {
        if (hs_coeff[i] >= 0) {
            minVal += hs_coeff[i] * mbr[i].first;
            maxVal += hs_coeff[i] * mbr[i].second;
        }
    }

    // Secondo ciclo: gestisci coefficienti negativi
    for (int i = 0; i < dims; ++i) {
        if (hs_coeff[i] < 0) {
            minVal += hs_coeff[i] * mbr[i].second;
            maxVal += hs_coeff[i] * mbr[i].first;
        }
    }

    if (maxVal < hs_known) {
        return BELOW;
    } else if (minVal > hs_known) {
        return ABOVE;
    } else {
        return OVERLAPPED;
    }
}

void QNode::appendHalfspace(const std::shared_ptr<HalfSpace>& hs) {
    PositionHS pos = MbrVersusHalfSpace(hs->coeff, hs->known);

    if (pos == BELOW) {
        covered.push_back(hs);
    } else if (pos == OVERLAPPED) {
        if (isLeaf()) {
            halfspaces.push_back(hs);

            if (halfspaces.size() > maxCapacity && level < maxLevel) {
                splitNode();
                for (const auto& oldHs : halfspaces) {
                    for (const auto& child : children) {
                        if (child->norm) {
                            child->appendHalfspace(oldHs);
                        }
                    }
                }
                halfspaces.clear();
            }
        } else {
            for (const auto& child : children) {
                if (child->norm) {
                    child->appendHalfspace(hs);
                }
            }
        }
    }
}

std::vector<std::shared_ptr<HalfSpace>> QNode::getTotalCovered() const {
    std::vector<std::shared_ptr<HalfSpace>> totalCovered = covered;
    QNode* ref = parent;

    while (ref != nullptr) {
        totalCovered.insert(totalCovered.end(), ref->covered.begin(), ref->covered.end());
        ref = ref->parent;
    }

    return totalCovered;
}

void QNode::splitNode() {
    ///if (level == maxLevel) return;
    if (!norm) return;

    auto subDivs = genSubdivisions();

    children.resize(numOfSubdivisions);
    for (int i = 0; i < numOfSubdivisions; ++i) {
        auto child = std::make_unique<QNode>(this, subDivs[i], dims);
        if (child->checkNodeValidity()) {
            child->norm = true;
        } else {
            child->norm = false;
        }
        children[i] = std::move(child);
    }

    leaf = false;
}

#define MAXDIM = 10

std::vector<std::vector<std::pair<double, double>>> QNode::genSubdivisions() {
    std::vector<double> mids(dims);
    for (int j = 0; j < dims; ++j) {
        mids[j] = (mbr[j].first + mbr[j].second) * 0.5;
    }

    std::vector<std::vector<std::pair<double, double>>> subdivisions(numOfSubdivisions);
    for (int i = 0; i < numOfSubdivisions; ++i) {
        subdivisions[i].resize(dims);
        for (int j = 0; j < dims; ++j) {
            if (i & (1 << j)) {
                subdivisions[i][j] = {mids[j], mbr[j].second};
            } else {
                subdivisions[i][j] = {mbr[j].first, mids[j]};
            }
        }
    }

    return subdivisions;
}

bool QNode::checkNodeValidity() {
    int numVertices = 1 << dims;
    for (int i = 0; i < numVertices; ++i) {
        double sum = 0.0;
        for (int j = 0; j < dims; ++j) {
            double coord = (i & (1 << j)) ? mbr[j].second : mbr[j].first;
            sum += coord;
        }
        if (sum <= normalizedMax) {
            return true;
        }
    }
    return false;
}