//
// Created by leona on 01/08/2024.
//

#ifndef QNODE_H
#define QNODE_H

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
    HalfSpace** covered;        // Halfspace coperti (array di puntatori)
    int numCovered;             // Numero di halfspace coperti
    HalfSpace** halfspaces;     // Halfspace nel nodo (array di puntatori)
    int numHalfspaces;          // Numero di halfspace nel nodo
    QNode** children;           // Puntatore ai figli
    int dims;                   // Numero di dimensioni del nodo

    // Funzioni membro
    inline bool isRoot() const { return parent == nullptr; }
    inline bool isLeaf() const { return leaf; }
    void setOrder();
    void insertHalfspaces(HalfSpace** hspaces, int numHSpaces);
    PositionHS MbrVersusHalfSpace(const double* hs_coeff, double hs_known);
    void clearHalfspaces();
    void splitNode();
    double*** genSubdivisions();
    bool checkNodeValidity();
    HalfSpace** getTotalCovered(int& totalCovered) const;
};

#endif // QNODE_H
