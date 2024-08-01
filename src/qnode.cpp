//
// Created by leona on 01/08/2024.
//

#include "qnode.h"
#include "halfspace.h"
#include <iostream>
#include <numeric>

#include "geom.h"

QNode::QNode(QNode* parent, const std::vector<std::array<double, 2>>& mbr)
    : mbr(mbr), norm(true), order(0), parent(parent) {}

bool QNode::isRoot() const {
    return parent == nullptr;
}

bool QNode::isLeaf() const {
    return children.empty();
}

int QNode::getOrder() {
    order = covered.size();
    QNode* ref = parent;

    while (ref && !ref->isRoot()) {
        order += ref->covered.size();
        ref = ref->parent;
    }

    return order;
}

std::vector<Halfspace> QNode::getCovered() {
    std::vector<Halfspace> coveredSpaces = covered;
    QNode* ref = parent;

    while (ref && !ref->isRoot()) {
        coveredSpaces.insert(coveredSpaces.end(), ref->covered.begin(), ref->covered.end());
        ref = ref->parent;
    }

    return coveredSpaces;
}

void QNode::insertHalfspaces(const std::vector<std::array<std::vector<double>, 2>>& masks, const std::vector<Halfspace>& halfspaces) {
    std::vector<double> incr(mbr.size());
    std::vector<double> half(mbr.size());

    for (size_t i = 0; i < mbr.size(); ++i) {
        incr[i] = (mbr[i][1] - mbr[i][0]) / 2.0;
        half[i] = (mbr[i][0] + mbr[i][1]) / 2.0;
    }

    const auto& ptsMask = masks[0];
    const auto& ndsMask = masks[1];

    std::vector<std::vector<double>> pts(ptsMask.size(), std::vector<double>(half.size()));

    for (size_t i = 0; i < ptsMask.size(); ++i) {
        for (size_t j = 0; j < half.size(); ++j) {
            pts[i][j] = incr[j] * ptsMask[i][j] + half[j];
        }
    }

    std::vector<std::vector<double>> coeff;
    std::vector<double> known;

    for (const auto& hs : halfspaces) {
        coeff.push_back(hs.coeff);
        known.push_back(hs.known);
    }

    std::vector<std::vector<int>> pos(pts.size(), std::vector<int>(halfspaces.size(), 0));

    for (size_t i = 0; i < pts.size(); ++i) {
        for (size_t j = 0; j < halfspaces.size(); ++j) {
            pos[i][j] = (std::inner_product(pts[i].begin(), pts[i].end(), coeff[j].begin(), 0.0) < known[j]) ? Position::IN : Position::OUT;
        }
    }

    for (size_t hs = 0; hs < halfspaces.size(); ++hs) {
        std::vector<int> rel;
        for (size_t i = 0; i < pos.size(); ++i) {
            if (pos[i][hs] != pos[0][hs]) {
                rel.push_back(i);
            }
        }

        std::vector<int> cross;
        for (size_t j = 0; j < ndsMask.size(); ++j) {
            bool crosses = false;
            for (const auto& r : rel) {
                if (ndsMask[r][j] > 0) {
                    crosses = true;
                    break;
                }
            }
            if (crosses) {
                cross.push_back(j);
            }
        }

        for (const auto& c : cross) {
            children[c]->halfspaces.push_back(halfspaces[hs]);
        }

        if (pos[0][hs] == Position::IN) {
            std::vector<int> notCross;
            for (size_t j = 0; j < ndsMask.size(); ++j) {
                bool crosses = false;
                for (const auto& r : rel) {
                    if (ndsMask[r][j] > 0) {
                        crosses = true;
                        break;
                    }
                }
                if (!crosses) {
                    notCross.push_back(j);
                }
            }

            for (const auto& nc : notCross) {
                children[nc]->covered.push_back(halfspaces[hs]);
            }
        }
    }
}