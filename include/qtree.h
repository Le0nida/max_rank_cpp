//
// Created by leona on 01/08/2024.
//

#ifndef QTREE_H
#define QTREE_H

#include "qnode.h"

class QTree {
public:
    QTree(int dims, int maxhsnode);

    void inserthalfspaces(const std::vector<long int>& halfspacesToInsert);
    QNode** getleaves(int& numLeaves);

private:
    int dims;        // Dimensionality of the space wrapped by the tree
    int maxhsnode;   // Maximum number of halfspaces a node can contain before being split up
    QNode* root;     // Reference to root node

    void createroot();
};

#endif // QTREE_H


