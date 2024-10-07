//
// Created by leona on 01/08/2024.
//

#include "qtree.h"
#include "qnode.h"
#include "geom.h"
#include "halfspace.h"
#include <vector>
#include <array>
#include <iostream>
#include <numeric>

extern int numOfSubdivisions;

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
    root->splitNode();
    return root;
}

// Insert halfspaces into the tree
void QTree::inserthalfspaces(const std::vector<long int>& halfspaces) {
    if (halfspaces.empty()) return;

    // Start with the root node
    root->insertHalfspaces(halfspaces);
}

// Retrieve all leaves of the QTree
std::vector<QNode*> QTree::getleaves() {
    std::vector<QNode*> leaves;
    std::vector<QNode*> to_search = {root};  // Contiene i nodi da elaborare

    while (!to_search.empty()) {
        QNode* current = to_search.back();
        to_search.pop_back();

        if (current->isNorm()) {
            if (current->isLeaf()) {
                leaves.push_back(current);  // Se è una foglia, aggiungila alla lista
            } else {
                for (int i = 0; i < numOfSubdivisions; i++) {
                    QNode* child = current->children[i];
                    // Verifica che il puntatore al figlio non sia nullo
                    if (!child) {
                        continue;
                    }
                    to_search.push_back(child);
                }
            }
        }
    }

    return leaves;
}
