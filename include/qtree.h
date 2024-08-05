//
// Created by leona on 01/08/2024.
//

#ifndef QTREE_H
#define QTREE_H

#include <vector>
#include <array>
#include <memory>
#include "qnode.h"
#include "geom.h"

class QTree {
public:
    QTree(int dims, int maxhsnode);
    ~QTree() = default;

    void insertHalfspaces(const std::vector<Halfspace>& halfspaces);
    std::vector<QNode*> getLeaves();

private:
    int dims;
    int maxhsnode;
    std::array<std::vector<std::vector<double>>, 2> masks;
    std::unique_ptr<QNode> root;

    std::unique_ptr<QNode> createRoot();
    void splitNode(QNode* node);
};

#endif // QTREE_H

