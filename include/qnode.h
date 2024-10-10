//
// Created by leona on 01/08/2024.
//

#ifndef QNODE_H
#define QNODE_H

#include <vector>

#include "halfspace.h"

extern int globalNodeID;
extern int numOfSubdivisions;
extern int normalizedMax;
extern int maxCapacity;

// Definizione dell'enumerazione per la posizione rispetto all'halfspace
enum PositionHS { BELOW, ABOVE, OVERLAPPED };

class QNode {
public:
    // Costruttore e distruttore
    QNode(QNode* parent, double** mbr, int dims);
    ~QNode();

    // Proprietà pubbliche
    long int nodeID;            // ID univoco del nodo
    QNode* parent;              // Puntatore al nodo genitore
    double** mbr;               // Minimum bounding region (array C-like)
    bool norm;                  // Flag di normalizzazione
    bool leaf;                  // Flag di foglia
    size_t order;               // Ordine del nodo
    std::vector<HalfSpace *> covered;
    std::vector<HalfSpace *> halfspaces;          // ID degli halfpsaces nel nodo
    QNode** children;           // Puntatore ai figli
    int dims;                   // Numero di dimensioni del nodo

    // Funzioni membro
    [[nodiscard]] inline bool isRoot() const { return parent == nullptr; }
    [[nodiscard]] inline bool isLeaf() const { return leaf; }
    void setOrder();
    PositionHS MbrVersusHalfSpace(const double* hs_coeff, double hs_known);
    void splitNode();
    double*** genSubdivisions();
    bool checkNodeValidity();
    std::vector<HalfSpace *> getTotalCovered(int& totalCovered) const;
    void appendHalfspace(HalfSpace * hsID);
};

#endif // QNODE_H
