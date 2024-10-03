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
#include <memory>

class Context; // Forward declaration

class QTree {
public:
    QTree(int dims, int maxhsnode, Context& ctx);
    void inserthalfspaces(const std::vector<long int>& halfspaces);
    std::vector<std::shared_ptr<QNode>> getleaves();

    bool hasHalfspaceBeenInserted(long int halfspaceID);
    void clearHalfspaceBeenInserted();

private:
    Context& ctx;
    int dims;
    int maxhsnode;
    std::array<std::vector<std::vector<double>>, 2> masks;
    std::shared_ptr<QNode> root;

    std::shared_ptr<QNode> createroot();

    std::unordered_set<long int> insertedHalfspaces;
};

#endif // QTREE_H

