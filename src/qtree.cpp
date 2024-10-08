//
// Created by leona on 01/08/2024.
//

#include "qtree.h"
#include <cstdlib> // For malloc, free

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

void QTree::inserthalfspaces(HalfSpace** halfspaces, int numHalfspaces) {
    if (numHalfspaces == 0) return;

    // Start with the root node
    root->insertHalfspaces(halfspaces, numHalfspaces);
}

QNode** QTree::getleaves(int& numLeaves) {
    QNode** leaves = nullptr;
    numLeaves = 0;
    QNode** to_search = (QNode**)malloc(sizeof(QNode*));
    int numToSearch = 1;
    to_search[0] = root;

    while (numToSearch > 0) {
        QNode* current = to_search[numToSearch - 1];
        numToSearch--;

        if (current->norm) {
            if (current->isLeaf()) {
                numLeaves++;
                leaves = (QNode**)realloc(leaves, numLeaves * sizeof(QNode*));
                leaves[numLeaves - 1] = current;
            } else {
                for (int i = 0; i < numOfSubdivisions; i++) {
                    QNode* child = current->children[i];
                    // Verify that the child pointer is not null
                    if (child) {
                        numToSearch++;
                        to_search = (QNode**)realloc(to_search, numToSearch * sizeof(QNode*));
                        to_search[numToSearch - 1] = child;
                    }
                }
            }
        }
    }
    free(to_search);

    return leaves;
}