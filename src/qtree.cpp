//
// Created by leona on 01/08/2024.
//

#include "qtree.h"
#include <vector>

extern int numOfSubdivisions;

QTree::QTree(int dims, int maxhsnode)
    : dims(dims), maxhsnode(maxhsnode) {
    createroot();
}

void QTree::createroot() {
    std::vector<std::pair<double, double>> mbr(dims, {0.0, 1.0});
    root = std::make_unique<QNode>(nullptr, mbr, dims);
    root->splitNode();
}

void QTree::inserthalfspaces(const std::vector<std::shared_ptr<HalfSpace>>& halfspacesToInsert) {
    for (const auto& hs : halfspacesToInsert) {
        for (auto& child : root->children) {
            if (child->norm) {
                child->appendHalfspace(hs);
            }
        }
    }
}

std::vector<QNode*> QTree::getleaves() {
    std::vector<QNode*> to_search = {root.get()};
    std::vector<QNode*> leaves;

    while (!to_search.empty()) {
        QNode* current = to_search.back();
        to_search.pop_back();

        if (current->norm) {
            if (current->isLeaf()) {
                leaves.push_back(current);
            } else {
                for (auto& child : current->children) {
                    if (child) {
                        to_search.push_back(child.get());
                    }
                }
            }
        }
    }

    return leaves;
}