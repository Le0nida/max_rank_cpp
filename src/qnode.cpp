//
// Created by leona on 01/08/2024.
//
#include "qnode.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstring>
#include <numeric>
#include "geom.h"
#include "qtree.h"

extern int numOfSubdivisions;
int normalizedMax = 1;
int maxCapacity = 10; // Capacità massima del nodo

// Contatore globale per assegnare ID unici
int globalNodeID = 0;

// Costruttore di QNode
QNode::QNode(QNode* parent, const std::vector<std::array<double, 2>>& mbr)
    : nodeID(globalNodeID++), parent(parent), mbr(mbr), norm(true), leaf(true), order(0), children(nullptr) {}

QNode::~QNode() {
    // Dealloca i nodi figli ricorsivamente
    if (children) {
        for (int i = 0; i < numOfSubdivisions; ++i) {
            if (children[i]) {
                delete children[i];
            }
        }
        delete[] children;
        children = nullptr;
    }
}

// Verifica se il nodo è la radice
bool QNode::isRoot() const {
    return parent == nullptr;
}

// Verifica se il nodo è una foglia
bool QNode::isLeaf() const {
    return leaf;
}

// Imposta l'ordine del nodo
void QNode::setOrder() {
    size_t localOrder = covered.size();
    QNode* ref = parent;

    while (ref != nullptr) {
        localOrder += ref->covered.size();
        ref = ref->parent;
    }

    order = localOrder;
}

// Ottiene gli halfspace coperti
std::vector<long int> QNode::getCovered() const {
    std::vector<long int> coveredSpaces = covered;
    QNode* ref = parent;

    while (ref != nullptr) {
        coveredSpaces.insert(coveredSpaces.end(), ref->covered.begin(), ref->covered.end());
        ref = ref->parent;
    }

    return coveredSpaces;
}

// Determina la posizione del MBR rispetto a un halfspace
PositionHS QNode::MbrVersusHalfSpace(const std::vector<double>& hs_coeff, double hs_known) {
    if (mbr.empty()) {
        // Gestisci il caso in cui mbr è vuoto
        return PositionHS::OVERLAPPED; // O qualsiasi comportamento predefinito tu preferisca
    }

    double minVal = 0.0;
    double maxVal = 0.0;
    size_t dims = hs_coeff.size();

    // Calcola i valori minimo e massimo dell'halfspace sul MBR
    for (size_t i = 0; i < dims; ++i) {
        double coeff = hs_coeff[i];
        if (coeff >= 0) {
            minVal += coeff * mbr[i][0]; // Usa il valore minimo del MBR per coefficienti positivi
            maxVal += coeff * mbr[i][1]; // Usa il valore massimo del MBR per coefficienti positivi
        } else {
            minVal += coeff * mbr[i][1]; // Usa il valore massimo del MBR per coefficienti negativi
            maxVal += coeff * mbr[i][0]; // Usa il valore minimo del MBR per coefficienti negativi
        }
    }

    // Confronta i valori min e max con hs_known per determinare la posizione
    if (maxVal < hs_known) {
        return PositionHS::BELOW;
    } else if (minVal > hs_known) {
        return PositionHS::ABOVE;
    } else {
        return PositionHS::OVERLAPPED;
    }
}

// Inserisce gli halfspace nel nodo e li distribuisce ai nodi figli se necessario
void QNode::insertHalfspaces(const std::vector<long int>& halfspaces) {
    size_t num_halfspaces = halfspaces.size();
    for (size_t i = 0; i < num_halfspaces; ++i) {
        long int hsID = halfspaces[i];
        auto hs = halfspaceCache->get(hsID);
        const std::vector<double>& hs_coeff = hs->coeff;
        double hs_known = hs->known;

        // Determina la posizione dell'halfspace rispetto al MBR del nodo
        PositionHS pos = MbrVersusHalfSpace(hs_coeff, hs_known);

        if (pos == PositionHS::BELOW) {
            // L'halfspace è completamente sotto il MBR, memorizzalo nei covered
            covered.push_back(hsID);
        } else if (pos == PositionHS::OVERLAPPED) {
            if (isLeaf()) {
                // Il nodo è una foglia, memorizza negli halfspaces
                this->halfspaces.push_back(hsID);

                // Verifica se è necessario suddividere il nodo
                if (this->halfspaces.size() > maxCapacity) {
                    // Suddividi il nodo
                    splitNode();

                    // Redistribuisci gli halfspaces tra i figli
                    for (int k = 0; k < numOfSubdivisions; k++) {
                        if (children[k] != nullptr) {
                            children[k]->insertHalfspaces(this->halfspaces);
                        }
                    }
                    // Svuota gli halfspaces da questo nodo
                    this->halfspaces.clear();
                }
            } else {
                // Il nodo non è una foglia, passa l'halfspace ai figli
                for (int k = 0; k < numOfSubdivisions; k++) {
                    if (children[k] != nullptr) {
                        children[k]->insertHalfspaces({hsID});
                    }
                }
            }
        }
    }
}

// Suddivide il nodo corrente in sotto-nodi
void QNode::splitNode() {
    if (!isNorm()) {
        // Non suddividere se il nodo è invalido
        return;
    }

    // Genera le suddivisioni (MBR dei figli)
    std::vector<std::vector<std::array<double, 2>>> subDivs = genSubdivisions();

    children = new QNode*[numOfSubdivisions];
    // Inizializza tutti i figli a nullptr
    for (int i = 0; i < numOfSubdivisions; ++i) {
        children[i] = nullptr;
    }

    int i = 0;
    for (const auto& child_mbr : subDivs) {
        // Crea un nuovo nodo figlio
        QNode* child = new QNode(this, child_mbr);

        // Verifica se il nodo figlio è valido
        if (child->checkNodeValidity()) {
            child->setNorm(true); // Imposta il nodo come valido
            children[i] = child;
        } else {
            // Il nodo è invalido, eliminato
            child->setNorm(false);
            delete child;
            children[i] = nullptr;
        }
        i++;
    }

    setLeaf(false); // Questo nodo non è più una foglia
}

// Genera le suddivisioni dell'MBR del nodo corrente
std::vector<std::vector<std::array<double, 2>>> QNode::genSubdivisions() {
    size_t dims = mbr.size();
    std::vector<std::vector<std::array<double, 2>>> subdivisions;

    size_t numOfSubdivisions = 1 << dims; // 2^dims combinazioni

    for (int i = 0; i < numOfSubdivisions; ++i) {
        std::vector<std::array<double, 2>> child_mbr(dims);

        for (size_t j = 0; j < dims; ++j) {
            double mid = (mbr[j][0] + mbr[j][1]) / 2.0; // Punto medio dell'MBR nella dimensione j

            if (i & (1 << j)) {
                // Usa la metà superiore nella dimensione j
                child_mbr[j][0] = mid;
                child_mbr[j][1] = mbr[j][1];
            } else {
                // Usa la metà inferiore nella dimensione j
                child_mbr[j][0] = mbr[j][0];
                child_mbr[j][1] = mid;
            }
        }

        subdivisions.push_back(child_mbr);
    }

    return subdivisions;
}

// Verifica se il nodo è valido in base al suo MBR
bool QNode::checkNodeValidity() {
    size_t dims = mbr.size();
    size_t num_vertices = 1 << dims; // Numero di vertici è 2^dims

    for (size_t i = 0; i < num_vertices; ++i) {
        double sum = 0.0;

        // Genera le coordinate del vertice
        for (size_t j = 0; j < dims; ++j) {
            double coord = (i & (1 << j)) ? mbr[j][1] : mbr[j][0]; // Usa il limite superiore o inferiore in base alla maschera
            sum += coord;
        }

        if (sum <= normalizedMax) {
            // Almeno un vertice è valido
            return true;
        }
    }
    // Tutti i vertici sono invalidi
    return false;
}

void QNode::clearHalfspaces() {
    halfspaces.clear();
}
