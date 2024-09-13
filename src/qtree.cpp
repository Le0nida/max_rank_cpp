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
    splitnode(root);
    globalCache.add(root);
    return root;
}

// Split the given node
void QTree::splitnode(std::shared_ptr<QNode> node) {
    const auto& mbr = node->getMBR();
    std::vector<double> mindim(dims);
    std::vector<double> maxdim(dims);

    for (int i = 0; i < dims; ++i) {
        mindim[i] = mbr[i][0];
        maxdim[i] = mbr[i][1];
    }

    for (int quad = 0; quad < (1 << dims); ++quad) {
        std::vector<double> child_mindim(dims);
        std::vector<double> child_maxdim(dims);

        for (int i = 0; i < dims; ++i) {
            if (quad & (1 << i)) {
                child_mindim[i] = (mindim[i] + maxdim[i]) / 2;
                child_maxdim[i] = maxdim[i];
            } else {
                child_mindim[i] = mindim[i];
                child_maxdim[i] = (mindim[i] + maxdim[i]) / 2;
            }
        }

        std::vector<std::array<double, 2>> child_mbr(dims);
        for (int i = 0; i < dims; ++i) {
            child_mbr[i] = {child_mindim[i], child_maxdim[i]};
        }

        auto child = std::make_unique<QNode>(node->getNodeID(), child_mbr);
        if (std::accumulate(child_mindim.begin(), child_mindim.end(), 0.0) >= 1.0) {
            child->setNorm(false);
        }
        node->addChildID(child->getNodeID());
        globalCache.add(std::move(child));
    }
}

// Insert halfspaces into the tree
void QTree::inserthalfspaces(const std::vector<long int>& halfspaces) {
    std::vector<int> to_search = {root->getNodeID()};  // Memorizza solo gli ID dei nodi
    root->setHalfspaces(halfspaces);

    while (!to_search.empty()) {
        int currentID = to_search.back();
        to_search.pop_back();

        std::shared_ptr<QNode> current = globalCache.get(currentID);  // Recupera il nodo dalla cache

        // Inserisce gli halfspaces e pulisce il nodo
        current->insertHalfspaces(masks, current->getHalfspaces());
        current->clearHalfspaces();

        // Copia locale sicura degli ID dei figli
        std::vector<long int> childrenIDs = current->getChildrenIDs();

        for (const auto childID : childrenIDs) {
            std::shared_ptr<QNode> child = globalCache.get(childID);

            if (child->isNorm()) {
                if (!child->isLeaf() && !child->getHalfspaces().empty()) {
                    to_search.push_back(childID);  // Aggiungi l'ID del figlio alla lista di ricerca
                } else if (child->getHalfspaces().size() > maxhsnode) {
                    splitnode(child);
                    to_search.push_back(childID);
                }
            }
        }
    }
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
