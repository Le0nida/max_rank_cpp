//
// Created by leona on 01/08/2024.
//

#include "qtree.h"
#include "qnode.h"
#include "lrucache.h"
#include "geom.h"
#include "halfspace.h"
#include <vector>
#include <array>
#include <iostream>
#include <numeric>

QTree::QTree(int dims, int maxhsnode, Context& ctx) : dims(dims), maxhsnode(maxhsnode), ctx(ctx) {
    auto masks_pair = genmasks(dims);
    masks = {masks_pair.first, masks_pair.second};
    root = createroot();
}

std::shared_ptr<QNode> QTree::createroot() {
    std::vector<std::array<double, 2>> mbr(dims, {0.0, 1.0});
    std::shared_ptr<QNode> root = std::make_shared<QNode>(ctx, -1, mbr);
    root->splitNode(ctx);
    ctx.cache->add(root);
    return root;
}

void QTree::inserthalfspaces(const std::vector<long int>& halfspaces) {
    if (halfspaces.empty()) return;

    root->insertHalfspaces(ctx, halfspaces);
}

std::vector<std::shared_ptr<QNode>> QTree::getleaves() {
    std::vector<std::shared_ptr<QNode>> leaves;
    std::vector<std::shared_ptr<QNode>> to_search = {root};

    while (!to_search.empty()) {
        std::shared_ptr<QNode> current = to_search.back();
        to_search.pop_back();

        ctx.cache->lockNode(current->getNodeID());

        if (current->isNorm()) {
            if (current->isLeaf()) {
                leaves.push_back(current);
            } else {
                std::vector<long int> childrenIDs = current->getChildrenIDs();
                for (const auto childID : childrenIDs) {
                    std::shared_ptr<QNode> child = ctx.cache->get(childID);
                    if (!child) {
                        std::cerr << "Error: Child with ID " << childID << " is null." << std::endl;
                        continue;
                    }
                    to_search.push_back(child);
                }
            }
        }

        ctx.cache->unlockNode(current->getNodeID());
    }

    return leaves;
}

bool QTree::hasHalfspaceBeenInserted(long int halfspaceID) {
    return insertedHalfspaces.find(halfspaceID) != insertedHalfspaces.end();
}

void QTree::clearHalfspaceBeenInserted() {
    insertedHalfspaces.clear();
}