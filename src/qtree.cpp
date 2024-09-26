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

#include "lrucache.h"

// Constructor for QTree
QTree::QTree(int dims, int maxhsnode) : dims(dims), maxhsnode(maxhsnode) {
    auto masks_pair = genmasks(dims);
    masks = {masks_pair.first, masks_pair.second};
    root = createroot();
}

// Create the root node and split it
std::shared_ptr<QNode> QTree::createroot() {
    std::vector<std::array<double, 2>> mbr(dims, {0.0, 1.0});
    std::shared_ptr<QNode> root = std::make_shared<QNode>(-1, mbr);
    root->splitNode();
    globalCache.add(root);
    return root;
}

// Insert halfspaces into the tree
void QTree::inserthalfspaces(const std::vector<long int>& halfspaces) {
    if (halfspaces.empty()) return;

    // Start with the root node
    root->insertHalfspaces(halfspaces);
}

// Retrieve all leaves of the QTree
std::vector<std::shared_ptr<QNode>> QTree::getleaves() {
    std::vector<std::shared_ptr<QNode>> leaves;
    std::vector<std::shared_ptr<QNode>> to_search = {root};  // Contiene i nodi da elaborare

    while (!to_search.empty()) {
        std::shared_ptr<QNode> current = to_search.back();
        to_search.pop_back();

        globalCache.lockNode(current->getNodeID());  // Blocca il nodo durante l'elaborazione

        if (current->isNorm()) {
            if (current->isLeaf()) {
                leaves.push_back(current);  // Se è una foglia, aggiungila alla lista
            } else {
                std::vector<long int> childrenIDs = current->getChildrenIDs();  // Copia locale sicura
                for (const auto childID : childrenIDs) {
                    std::shared_ptr<QNode> child = globalCache.get(childID);

                    // Verifica che il puntatore al figlio non sia nullo
                    if (!child) {
                        std::cerr << "Error: Child with ID " << childID << " is null." << std::endl;
                        continue;
                    }
                    to_search.push_back(child);
                }
            }
        }

        globalCache.unlockNode(current->getNodeID());  // Sblocca il nodo dopo l'elaborazione
    }

    return leaves;
}
