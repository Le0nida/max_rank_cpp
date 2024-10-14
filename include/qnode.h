//
// Created by leona on 01/08/2024.
//

#ifndef QNODE_H
#define QNODE_H

#include <vector>
#include <memory>
#include "halfspace.h"

// Definizione dell'enumerazione per la posizione rispetto all'halfspace
enum PositionHS { BELOW, ABOVE, OVERLAPPED };

class QNode {
public:
    long int nodeID;
    QNode* parent;
    std::vector<std::pair<double, double>> mbr; // MBR rappresentato come vettore di coppie (min, max)
    bool norm;
    bool leaf;
    size_t order;
    std::vector<std::shared_ptr<HalfSpace>> covered;
    std::vector<std::shared_ptr<HalfSpace>> halfspaces;
    std::vector<std::unique_ptr<QNode>> children;
    int dims;

    QNode(QNode* parent, const std::vector<std::pair<double, double>>& mbr, int dims);

    void setOrder();
    PositionHS MbrVersusHalfSpace(const std::vector<double>& hs_coeff, double hs_known);
    void splitNode();
    std::vector<std::vector<std::pair<double, double>>> genSubdivisions();
    bool checkNodeValidity();
    [[nodiscard]] std::vector<std::shared_ptr<HalfSpace>> getTotalCovered() const;
    void appendHalfspace(const std::shared_ptr<HalfSpace>& hs);

    [[nodiscard]] bool isRoot() const { return parent == nullptr; }
    [[nodiscard]] bool isLeaf() const { return leaf; }
};


#endif // QNODE_H
