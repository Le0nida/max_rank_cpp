//
// Created by leona on 01/08/2024.
//

#include "qtree.h"
#include "qnode.h"
#include "halfspace.h"
#include "geom.h"
#include <vector>
#include <array>
#include <numeric>

// Constructor for QTree
QTree::QTree(int dims, int maxhsnode) : dims(dims), maxhsnode(maxhsnode) {
    auto masks_pair = genmasks(dims);
    masks = {masks_pair.first, masks_pair.second};
    root = createroot();
}

// Create the root node and split it
QNode* QTree::createroot() {
    std::vector<std::array<double, 2>> mbr(dims, {0.0, 1.0});
    auto* root = new QNode(nullptr, mbr);
    splitnode(root);
    return root;
}

// Split the given node
void QTree::splitnode(QNode* node) {
    const auto& mbr = node->getMBR();
    std::vector<double> mindim(dims);
    std::vector<double> maxdim(dims);

    for (int i = 0; i < dims; ++i) {
        mindim[i] = mbr[i][0];
        maxdim[i] = mbr[i][1];
    }

    for (int quad = 0; quad < (1 << dims); ++quad) {
        std::vector<double> child_mindim(dims);
        std::vector<double> child_maxdim(dims);

        for (int i = 0; i < dims; ++i) {
            if (quad & (1 << i)) {
                child_mindim[i] = (mindim[i] + maxdim[i]) / 2;
                child_maxdim[i] = maxdim[i];
            } else {
                child_mindim[i] = mindim[i];
                child_maxdim[i] = (mindim[i] + maxdim[i]) / 2;
            }
        }

        std::vector<std::array<double, 2>> child_mbr(dims);
        for (int i = 0; i < dims; ++i) {
            child_mbr[i] = {child_mindim[i], child_maxdim[i]};
        }

        auto child = std::make_unique<QNode>(node, child_mbr);
        if (std::accumulate(child_mindim.begin(), child_mindim.end(), 0.0) >= 1.0) {
            child->setNorm(false);
        }
        node->addChild(std::move(child));
    }
}

// Insert halfspaces into the tree
void QTree::inserthalfspaces(const std::vector<Halfspace>& halfspaces) {
    std::vector<QNode*> to_search = {root};
    root->setHalfspaces(halfspaces);

    while (!to_search.empty()) {
        QNode* current = to_search.back();
        to_search.pop_back();

        current->insertHalfspaces(masks, current->getHalfspaces());
        current->clearHalfspaces();

        for (const auto& child : current->getChildren()) {
            if (child->isNorm()) {
                if (!child->isLeaf() && !child->getHalfspaces().empty()) {
                    to_search.push_back(child.get());
                } else if (child->getHalfspaces().size() > maxhsnode) {
                    splitnode(child.get());
                    to_search.push_back(child.get());
                }
            }
        }
    }
}

// Retrieve all leaves of the QTree
std::vector<QNode*> QTree::getleaves() {
    std::vector<QNode*> leaves;
    std::vector<QNode*> to_search = {root};

    while (!to_search.empty()) {
        QNode* current = to_search.back();
        to_search.pop_back();

        if (current->isNorm()) {
            if (current->isLeaf()) {
                leaves.push_back(current);
            } else {
                for (const auto& child : current->getChildren()) {
                    to_search.push_back(child.get());
                }
            }
        }
    }

    return leaves;
}
