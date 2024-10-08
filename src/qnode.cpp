//
// Created by leona on 01/08/2024.
//

#include "qnode.h"
#include <cmath>
#include <iostream>
#include <cstring> // Per memcpy
#include <cstdlib> // Per malloc e free

// Variabili esterne
int normalizedMax = 1;
int maxCapacity = 10; // Capacità massima del nodo
int globalNodeID = 0;

QNode::QNode(QNode* parent, double** mbr, int dims)
    : parent(parent), mbr(mbr), norm(true), leaf(true), order(0),
      covered(nullptr), numCovered(0), halfspaces(nullptr), numHalfspaces(0),
      children(nullptr), dims(dims)
{
    nodeID = globalNodeID++;
}

QNode::~QNode() {
    // Dealloca gli halfspace coperti
    if (covered) {
        free(covered);
        covered = nullptr;
    }

    // Dealloca gli halfspace del nodo
    if (halfspaces) {
        free(halfspaces);
        halfspaces = nullptr;
    }

    // Dealloca i nodi figli ricorsivamente
    if (children) {
        for (int i = 0; i < numOfSubdivisions; ++i) {
            if (children[i]) {
                delete children[i];
            }
        }
        free(children);
        children = nullptr;
    }

    // Dealloca il mbr
    if (mbr) {
        for (int i = 0; i < dims; ++i) {
            free(mbr[i]);
        }
        free(mbr);
        mbr = nullptr;
    }
}

void QNode::setOrder() {
    size_t localOrder = numCovered;
    QNode* ref = parent;

    while (ref != nullptr) {
        localOrder += ref->numCovered;
        ref = ref->parent;
    }

    order = localOrder;
}

PositionHS QNode::MbrVersusHalfSpace(const double* hs_coeff, double hs_known) {
    double minVal = 0.0;
    double maxVal = 0.0;

    // Calcola i valori minimo e massimo dell'halfspace sul MBR
    for (int i = 0; i < dims; ++i) {
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
        return BELOW;
    } else if (minVal > hs_known) {
        return ABOVE;
    } else {
        return OVERLAPPED;
    }
}

void QNode::insertHalfspaces(HalfSpace** hspaces, int numHSpaces) {
    for (int i = 0; i < numHSpaces; ++i) {
        HalfSpace* hs = hspaces[i];
        const double* hs_coeff = hs->coeff;
        double hs_known = hs->known;

        // Determina la posizione dell'halfspace rispetto al MBR del nodo
        PositionHS pos = MbrVersusHalfSpace(hs_coeff, hs_known);

        if (pos == BELOW) {
            // L'halfspace è completamente sotto il MBR, memorizzalo nei covered
            numCovered++;
            covered = (HalfSpace**)realloc(covered, numCovered * sizeof(HalfSpace*));
            covered[numCovered - 1] = hs;
        } else if (pos == OVERLAPPED) {
            if (isLeaf()) {
                // Il nodo è una foglia, memorizza negli halfspaces
                numHalfspaces++;
                halfspaces = (HalfSpace**)realloc(halfspaces, numHalfspaces * sizeof(HalfSpace*));
                halfspaces[numHalfspaces - 1] = hs;

                // Verifica se è necessario suddividere il nodo
                if (numHalfspaces > maxCapacity) {
                    // Suddividi il nodo
                    splitNode();

                    // Redistribuisci gli halfspaces tra i figli
                    for (int k = 0; k < numOfSubdivisions; k++) {
                        if (children[k] != nullptr) {
                            children[k]->insertHalfspaces(halfspaces, numHalfspaces);
                        }
                    }
                    // Svuota gli halfspaces da questo nodo
                    clearHalfspaces();
                }
            } else {
                // Il nodo non è una foglia, passa l'halfspace ai figli
                for (int k = 0; k < numOfSubdivisions; k++) {
                    if (children[k] != nullptr) {
                        children[k]->insertHalfspaces(&hs, 1);
                    }
                }
            }
        }
    }
}

void QNode::clearHalfspaces() {
    if (halfspaces) {
        free(halfspaces);
        halfspaces = nullptr;
    }
    numHalfspaces = 0;
}

void QNode::splitNode() {
    if (!norm) {
        // Do not split if the node is invalid
        return;
    }

    int numSubs;
    double*** subDivs = genSubdivisions(numSubs);

    children = (QNode**)malloc(numSubs * sizeof(QNode*));
    // Initialize all children to nullptr
    for (int i = 0; i < numSubs; ++i) {
        children[i] = nullptr;
    }

    for (int i = 0; i < numSubs; ++i) {
        // Create a new child node
        QNode* child = new QNode(this, subDivs[i], dims);

        // Check if the child node is valid
        if (child->checkNodeValidity()) {
            child->norm = true; // Set the node as valid
            children[i] = child;
            subDivs[i] = nullptr; // The child now owns the mbr
        } else {
            // The node is invalid, delete it
            child->norm = false;
            delete child;
            children[i] = nullptr;
            subDivs[i] = nullptr; // The child has freed the mbr
        }
    }

    // Deallocate the subdivisions
    for (int i = 0; i < numSubs; ++i) {
        if (subDivs[i]) {
            for (int j = 0; j < dims; ++j) {
                if (subDivs[i][j]) {
                    free(subDivs[i][j]);
                }
            }
            free(subDivs[i]);
        }
    }
    free(subDivs);

    leaf = false; // This node is no longer a leaf
}


double*** QNode::genSubdivisions(int& numSubs) {
    numSubs = 1 << dims; // 2^dims combinazioni
    double*** subdivisions = (double***)malloc(numSubs * sizeof(double**));

    for (int i = 0; i < numSubs; ++i) {
        double** child_mbr = (double**)malloc(dims * sizeof(double*));
        for (int j = 0; j < dims; ++j) {
            child_mbr[j] = (double*)malloc(2 * sizeof(double));
            double mid = (mbr[j][0] + mbr[j][1]) * 0.5; // Punto medio dell'MBR nella dimensione j

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
        subdivisions[i] = child_mbr;
    }

    return subdivisions;
}

bool QNode::checkNodeValidity() {
    int num_vertices = 1 << dims; // Numero di vertici è 2^dims

    for (int i = 0; i < num_vertices; ++i) {
        double sum = 0.0;

        // Genera le coordinate del vertice
        for (int j = 0; j < dims; ++j) {
            double coord = (i & (1 << j)) ? mbr[j][1] : mbr[j][0]; // Usa il limite superiore o inferiore
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
