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
    std::vector<QNode*> getleaves();

private:
    int dims;  // Dimensionality of the space wrapped by the tree
    int maxhsnode;  // Maximum number of halfspaces a node can contain before being split up
    std::array<std::vector<std::vector<double>>, 2> masks;  // Masks used in halfspace insertion
    QNode* root;  // Reference to root node

    QNode* createroot();
};

#endif // QTREE_H

