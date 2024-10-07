//
// Created by leona on 01/08/2024.
//

#ifndef QNODE_H
#define QNODE_H

#include <vector>
#include <array>
#include "halfspace.h"

extern int globalNodeID;

enum class PositionHS { BELOW, ABOVE, OVERLAPPED };

class QNode {
public:
    QNode(QNode* parent, const std::vector<std::array<double, 2>>& mbr);
    ~QNode();

    // Disabilita il costruttore di copia e l'assegnazione tramite copia
    QNode(const QNode&) = delete;
    QNode& operator=(const QNode&) = delete;

    // Metodo per ottenere l'ID del nodo
    int getNodeID() const { return nodeID; }

    // Verifica se il nodo è la radice
    bool isRoot() const;

    // Verifica se il nodo è una foglia
    bool isLeaf() const;

    // Imposta l'ordine del nodo
    void setOrder();

    // Ottiene gli halfspace coperti
    std::vector<long int> getCovered() const;

    // Inserisce gli halfspace nel nodo
    void insertHalfspaces(const std::vector<long int>& halfspaces);

    // Determina la posizione del MBR rispetto a un halfspace
    PositionHS MbrVersusHalfSpace(const std::vector<double>& hs_coeff, double hs_known);

    // Getter e setter
    const std::vector<std::array<double, 2>>& getMBR() const { return mbr; }
    bool isNorm() const { return norm; }
    void setNorm(bool norm) { this->norm = norm; }
    void setLeaf(bool leaf) { this->leaf = leaf; }
    size_t getOrder() const { return order; }

    const std::vector<long int>& getHalfspaces() const { return halfspaces; }
    void clearHalfspaces();

    void splitNode();
    std::vector<std::vector<std::array<double, 2>>> genSubdivisions();
    bool checkNodeValidity();

    QNode** children;  // Puntatore ai figli

private:
    long int nodeID;                        // ID univoco del nodo
    QNode* parent;                          // Puntatore al nodo genitore

    std::vector<std::array<double, 2>> mbr; // Minimum bounding region
    bool norm;                              // Flag di normalizzazione
    bool leaf;                              // Flag di foglia
    size_t order;                           // Ordine del nodo

    std::vector<long int> covered;          // Halfspace coperti
    std::vector<long int> halfspaces;       // Halfspace nel nodo
};

#endif // QNODE_H
