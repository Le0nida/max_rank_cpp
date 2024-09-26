//
// Created by leona on 01/08/2024.
//

#ifndef QTREE_H
#define QTREE_H

#include "qnode.h"
#include "halfspace.h"
#include <vector>
#include <array>
#include <unordered_set>

class QTree {
public:
    QTree(int dims, int maxhsnode);
    void inserthalfspaces(const std::vector<long int>& halfspaces);
    std::vector<std::shared_ptr<QNode>> getleaves();

    bool hasHalfspaceBeenInserted(long int halfspaceID) {
        return insertedHalfspaces.find(halfspaceID) != insertedHalfspaces.end();
    }
    void clearHalfspaceBeenInserted() {
        insertedHalfspaces.clear();
    }

private:
    int dims;  // Dimensionality of the space wrapped by the tree
    int maxhsnode;  // Maximum number of halfspaces a node can contain before being split up
    std::array<std::vector<std::vector<double>>, 2> masks;  // Masks used in halfspace insertion
    std::shared_ptr<QNode> root;  // Reference to root node

    std::shared_ptr<QNode> createroot();

    std::unordered_set<long int> insertedHalfspaces;
};

#endif // QTREE_H

