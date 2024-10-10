//
// Created by leona on 01/08/2024.
//

#include "qtree.h"
#include <cstdlib> // For malloc, free
#include <string.h>
#include <vector>

extern int numOfSubdivisions;

// Constructor for QTree
QTree::QTree(int dims, int maxhsnode) : dims(dims), maxhsnode(maxhsnode) {
    createroot();
}

void QTree::createroot() {
    // Dynamic allocation for MBR (C-style 2D array)
    double** mbr = (double**)malloc(dims * sizeof(double*));

    for (int i = 0; i < dims; ++i) {
        mbr[i] = (double*)malloc(2 * sizeof(double));
        mbr[i][0] = 0.0;  // Lower bound
        mbr[i][1] = 1.0;  // Upper bound
    }

    // Create root with dynamically allocated MBR
    root = new QNode(nullptr, mbr, dims);

    // Split the root node
    root->splitNode();
}

// IN input passso la lista di id da inserire
void QTree::inserthalfspaces(const std::vector<HalfSpace *>& halfspacesToInsert) {
    if (halfspacesToInsert.empty()) return;

    for (auto hsID: halfspacesToInsert)
    {
        // Start with the root node
        for (int j = 0; j < numOfSubdivisions; j++) {                     //to each of the 2^d children of root
            if (root->children[j]->norm)
                root->children[j]->appendHalfspace(hsID);
        }
    }
}

QNode** QTree::getleaves(int& numLeaves) {
    std::vector<QNode*> to_search;
    std::vector<QNode*> leaves;

    to_search.push_back(root);

    while (!to_search.empty()) {
        QNode* current = to_search.back();
        to_search.pop_back();

        if (current->norm) {
            if (current->isLeaf()) {
                leaves.push_back(current);
            } else {
                for (int i = 0; i < numOfSubdivisions; i++) {
                    QNode* child = current->children[i];
                    if (child) {
                        to_search.push_back(child);
                    }
                }
            }
        }
    }

    // Copy leaves from std::vector to a dynamically allocated array
    numLeaves = static_cast<int>(leaves.size());
    QNode** result = (QNode**)malloc(numLeaves * sizeof(QNode*));
    memcpy(result, leaves.data(), numLeaves * sizeof(QNode*));

    return result;
}