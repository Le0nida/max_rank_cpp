//
// Created by leona on 01/08/2024.
//

#include "qnode.h"
#include "halfspace.h"
#include <iostream>
#include <fstream>
#include <cstring>  // Per memcpy
#include "lrucache.h"
#include <numeric>
#include "geom.h"

// Counter globale per assegnare ID unici (opzione semplice)
int globalNodeID = 0;

// Constructor for QNode, initializes the member variables.
QNode::QNode(int parentID, const std::vector<std::array<double, 2>>& mbr)
    : mbr(mbr), norm(true), order(0), parentID(parentID), nodeID(globalNodeID++){}

// Checks if the node is the root of the tree.
bool QNode::isRoot() const {
    return parentID == -1;
}

// Checks if the node is a leaf (i.e., has no children).
bool QNode::isLeaf() const {
    return childrenIDs.empty();
}

// Computes the order of the node by traversing back up the tree.
void QNode::setOrder() const {
    size_t localOrder = covered.size();  // Initialize order with the size of covered halfspaces in the current node.
    std::shared_ptr<QNode> ref = globalCache.get(parentID);

    // Traverse back up the tree to accumulate the order from parent nodes.
    while (ref->parentID >= 0) {
        localOrder += ref->covered.size();
        ref = globalCache.get(ref->parentID);
    }

    order = localOrder;
}

void QNode::clearHalfspaces()
{
    if (!halfspaces.empty()) {
        halfspaces.clear();
    } else {
        std::cerr << "Warning: Attempting to clear empty halfspaces." << std::endl;
    }
}


// Retrieves the covering halfspaces by traversing back up the tree.
std::vector<long int> QNode::getCovered() const{
    std::vector<long int> coveredSpaces = covered;  // Start with halfspaces covered in the current node.
    std::shared_ptr<QNode> ref = globalCache.get(parentID);;

    // Traverse back up the tree to accumulate the covered halfspaces from parent nodes.
    while (ref && !ref->isRoot()) {
        coveredSpaces.insert(coveredSpaces.end(), ref->covered.begin(), ref->covered.end());
        ref = globalCache.get(ref->parentID);;
    }

    return coveredSpaces;
}

// Inserts halfspaces into the node and distributes them to children nodes if necessary.
void QNode::insertHalfspaces(const std::array<std::vector<std::vector<double>>, 2>& masks, const std::vector<long int>& halfspaces) {
    const size_t mbr_size = mbr.size();
    double incr[mbr_size];
    double half[mbr_size];

    // Precalcolo di incr e half
    for (size_t i = 0; i < mbr_size; ++i) {
        incr[i] = (mbr[i][1] - mbr[i][0]) / 2.0;
        half[i] = (mbr[i][0] + mbr[i][1]) / 2.0;
    }

    if (masks.empty()) {
        std::cout << "Masks are empty, exiting function." << std::endl;
        return;
    }
    const auto& ptsMask = masks[0];  // Point masks for halfspace insertion
    const auto& ndsMask = masks[1];  // Node masks for halfspace insertion

    const size_t num_pts = ptsMask.size();
    const size_t num_halfspaces = halfspaces.size();

    // Calcola i punti per le maschere
    double pts[num_pts][mbr_size];
    for (size_t i = 0; i < num_pts; ++i) {
        for (size_t j = 0; j < mbr_size; ++j) {
            pts[i][j] = incr[j] * ptsMask[i][j] + half[j];
        }
    }

    // Estrai coefficienti e valori noti dalle halfspaces in array statici
    double coeff[num_halfspaces][mbr_size];
    double known[num_halfspaces];
    for (size_t i = 0; i < num_halfspaces; ++i) {
        auto hs = halfspaceCache->get(halfspaces[i]);
        memcpy(coeff[i], hs->coeff.data(), mbr_size * sizeof(double));
        known[i] = hs->known;
    }

    // Determina la posizione dei punti relativi a ciascuna halfspace
    int pos[num_pts][num_halfspaces];
    for (size_t i = 0; i < num_pts; ++i) {
        for (size_t j = 0; j < num_halfspaces; ++j) {
            double dot_product = 0.0;
            for (size_t k = 0; k < mbr_size; ++k) {
                dot_product += pts[i][k] * coeff[j][k];
            }
            pos[i][j] = (dot_product < known[j]) ? static_cast<int>(Position::IN) : static_cast<int>(Position::OUT);
        }
    }

    // Distribuire le halfspaces ai nodi figli in base alle loro posizioni
    for (size_t hs = 0; hs < num_halfspaces; ++hs) {
        std::vector<size_t> rel;  // Posizioni relative dove le posizioni dei punti differiscono
        for (size_t i = 1; i < num_pts; ++i) {  // Partendo da 1, evitiamo il confronto con sé stessi
            if (pos[i][hs] != pos[0][hs]) {
                rel.push_back(i);
            }
        }

        // Determina quali nodi figli la halfspace attraversa
        std::vector<size_t> cross;
        for (size_t j = 0; j < ndsMask[0].size(); ++j) {
            bool crosses = false;
            for (const auto& r : rel) {
                if (ndsMask[r][j] > 0) {
                    crosses = true;
                    break;
                }
            }
            if (crosses) {
                cross.push_back(j);
            }
        }

        // Aggiungi halfspace ai nodi figli appropriati
        for (const auto& c : cross) {
            if (c < childrenIDs.size()) {
                int childID = childrenIDs[c];
                auto child = globalCache.get(childID);  // Recupera il figlio dalla cache usando l'ID
                child->halfspaces.push_back(halfspaces[hs]);  // Inserisci l'halfspace nel nodo figlio
            } else {
                std::cerr << "Index out of bounds: " << c << " >= " << childrenIDs.size() << std::endl;
                return;
            }
        }

        // Aggiungi halfspace alla lista covered dei nodi figli appropriati se non attraversa
        if (pos[0][hs] == static_cast<int>(Position::IN)) {
            int sum_mask[ndsMask[0].size()] = {0};

            // Somma i valori in nds_mask[rel]
            for (size_t r : rel) {
                for (size_t j = 0; j < ndsMask[r].size(); ++j) {
                    sum_mask[j] += ndsMask[r][j];
                }
            }

            // Trova gli indici dove la somma è zero
            for (size_t nc = 0; nc < ndsMask[0].size(); ++nc) {
                if (sum_mask[nc] == 0) {
                    if (nc < childrenIDs.size()) {
                        int childID = childrenIDs[nc];
                        std::shared_ptr<QNode> child = globalCache.get(childID);  // Recupera il figlio dalla cache usando l'ID
                        child->covered.push_back(halfspaces[hs]);  // Aggiungi alla lista covered
                    } else {
                        std::cerr << "Index out of bounds: " << nc << " >= " << childrenIDs.size() << std::endl;
                        return;
                    }
                }
            }
        }
    }
}


// Serializza l'oggetto QNode su disco
void QNode::saveToDisk(const std::string& filePath) {
    std::ofstream out(filePath, std::ios::binary);
    if (out.is_open()) {
        // Serialize the nodeID and parentID
        out.write(reinterpret_cast<const char*>(&nodeID), sizeof(nodeID));
        out.write(reinterpret_cast<const char*>(&parentID), sizeof(parentID));

        // Serialize the MBR
        size_t mbrSize = mbr.size();
        out.write(reinterpret_cast<const char*>(&mbrSize), sizeof(mbrSize));
        out.write(reinterpret_cast<const char*>(mbr.data()), mbrSize * sizeof(std::array<double, 2>));

        // Serialize the norm flag
        out.write(reinterpret_cast<const char*>(&norm), sizeof(norm));

        // Serialize the children IDs
        size_t childrenCount = childrenIDs.size();
        out.write(reinterpret_cast<const char*>(&childrenCount), sizeof(childrenCount));
        out.write(reinterpret_cast<const char*>(childrenIDs.data()), childrenCount * sizeof(long int));

        // Serialize the halfspaces
        size_t halfspaceCount = halfspaces.size();
        out.write(reinterpret_cast<const char*>(&halfspaceCount), sizeof(halfspaceCount));
        out.write(reinterpret_cast<const char*>(halfspaces.data()), halfspaceCount * sizeof(long int));

        // Serialize the covered halfspaces
        size_t coveredCount = covered.size();
        out.write(reinterpret_cast<const char*>(&coveredCount), sizeof(coveredCount));
        out.write(reinterpret_cast<const char*>(covered.data()), coveredCount * sizeof(long int));

        out.close();
    } else {
        std::cerr << "Failed to open file for writing: " << filePath << std::endl;
    }
}


// Carica l'oggetto QNode da disco
void QNode::loadFromDisk(const std::string& filePath) {
    std::ifstream in(filePath, std::ios::binary);
    if (in.is_open()) {
        // Deserialize the nodeID and parentID
        in.read(reinterpret_cast<char*>(&nodeID), sizeof(nodeID));
        in.read(reinterpret_cast<char*>(&parentID), sizeof(parentID));

        // Deserialize the MBR
        size_t mbrSize;
        in.read(reinterpret_cast<char*>(&mbrSize), sizeof(mbrSize));
        mbr.resize(mbrSize);
        in.read(reinterpret_cast<char*>(mbr.data()), mbrSize * sizeof(std::array<double, 2>));

        // Deserialize the norm flag
        in.read(reinterpret_cast<char*>(&norm), sizeof(norm));

        // Deserialize the children IDs
        size_t childrenCount;
        in.read(reinterpret_cast<char*>(&childrenCount), sizeof(childrenCount));
        childrenIDs.resize(childrenCount);
        in.read(reinterpret_cast<char*>(childrenIDs.data()), childrenCount * sizeof(long int));

        // Deserialize the halfspaces
        size_t halfspaceCount;
        in.read(reinterpret_cast<char*>(&halfspaceCount), sizeof(halfspaceCount));
        halfspaces.resize(halfspaceCount);
        in.read(reinterpret_cast<char*>(halfspaces.data()), halfspaceCount * sizeof(long int));

        // Deserialize the covered halfspaces
        size_t coveredCount;
        in.read(reinterpret_cast<char*>(&coveredCount), sizeof(coveredCount));
        covered.resize(coveredCount);
        in.read(reinterpret_cast<char*>(covered.data()), coveredCount * sizeof(long int));

        in.close();
    } else {
        std::cerr << "Failed to open file for reading: " << filePath << std::endl;
    }
}
