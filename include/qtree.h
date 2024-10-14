//
// Created by leona on 01/08/2024.
//

#ifndef QTREE_H
#define QTREE_H

#include "qnode.h"
#include <memory>

class QTree {
public:
    QTree(int dims, int maxhsnode);

    void inserthalfspaces(const std::vector<std::shared_ptr<HalfSpace>>& halfspacesToInsert);
    std::vector<QNode*> getleaves();

private:
    int dims;
    int maxhsnode;
    std::unique_ptr<QNode> root;

    void createroot();
};

#endif // QTREE_H
