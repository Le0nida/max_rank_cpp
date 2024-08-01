//
// Created by leona on 01/08/2024.
//

#include "qtree.h"
#include "halfspace.h"
#include <iostream>
#include <numeric>
#include <bitset>

QTree::QTree(int dims, int maxhsnode)
    : dims(dims), maxhsnode(maxhsnode), masks(genMasks(dims)) {
    root = createRoot();
}

std::unique_ptr<QNode> QTree::createRoot() {
    std::vector<std::array<double, 2>> mbr(dims, {0.0, 1.0});
    auto root = std::make_unique<QNode>(nullptr, mbr);
    splitNode(root.get());
    return root;
}

void QTree::splitNode(QNode* node) {
    const auto& mbr = node->getMBR();
    const auto& mindim = mbr[0];
    const auto& maxdim = mbr[1];

    for (int quad = 0; quad < (1 << dims); ++quad) {
        std::bitset<32> qbin(quad); // Assuming dims <= 32

        std::array<double, 2> child_mindim;
        std::array<double, 2> child_maxdim;

        for (size_t i = 0; i < mindim.size(); ++i) {
            if (qbin[i] == 0) {
                child_mindim[i] = mindim[i];
                child_maxdim[i] = (mindim[i] + maxdim[i]) / 2;
            } else {
                child_mindim[i] = (mindim[i] + maxdim[i]) / 2;
                child_maxdim[i] = maxdim[i];
            }
        }

        std::vector<std::array<double, 2>> child_mbr = {child_mindim, child_maxdim};
        auto child = std::make_unique<QNode>(node, child_mbr);

        if (std::accumulate(child_mindim.begin(), child_mindim.end(), 0.0) >= 1) {
            child->setNorm(false);
        }

        node->addChild(std::move(child));
    }
}

void QTree::insertHalfspaces(const std::vector<Halfspace>& halfspaces) {
    std::vector<QNode*> to_search = {root.get()};
    root->setHalfspaces(halfspaces);

    while (!to_search.empty()) {
        QNode* current = to_search.back();
        to_search.pop_back();

        current->insertHalfspaces(masks, current->getHalfspaces());
        current->clearHalfspaces();

        for (auto& child : current->getChildren()) {
            if (child->isNorm()) {
                if (!child->isLeaf() && !child->getHalfspaces().empty()) {
                    to_search.push_back(child.get());
                } else if (child->getHalfspaces().size() > maxhsnode) {
                    splitNode(child.get());
                    to_search.push_back(child.get());
                }
            }
        }
    }
}

std::vector<QNode*> QTree::getLeaves() {
    std::vector<QNode*> leaves;
    std::vector<QNode*> to_search = {root.get()};

    while (!to_search.empty()) {
        QNode* current = to_search.back();
        to_search.pop_back();

        if (current->isNorm()) {
            if (current->isLeaf()) {
                leaves.push_back(current);
            } else {
                for (auto& child : current->getChildren()) {
                    to_search.push_back(child.get());
                }
            }
        }
    }

    return leaves;
}