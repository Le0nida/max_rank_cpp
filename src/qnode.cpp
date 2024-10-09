//
// Created by leona on 01/08/2024.
//

#include "qnode.h"
#include <cmath>
#include <iostream>
#include <cstring> // Per memcpy
#include <cstdlib> // Per malloc e free
#include <map>

// Variabili esterne
int normalizedMax = 1;
int maxCapacity = 10; // Capacità massima del nodo
int globalNodeID = 0;
extern std::map<long int, HalfSpace*> HalfSpacesMap;;

QNode::QNode(QNode* parent, double** mbr, int dims)
    : parent(parent), mbr(mbr), norm(true), leaf(true), order(0),
      covered(nullptr), numCovered(0), halfspaces({}),
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

void QNode::appendHalfspace(long int hsID)
{
    PositionHS pos = MbrVersusHalfSpace(HalfSpacesMap[hsID]->coeff, HalfSpacesMap[hsID]->known);
    // Determina la posizione dell'halfspace rispetto al MBR del nodo
    if (pos == BELOW) {
        // L'halfspace è completamente sotto il MBR, memorizzalo nei covered
        numCovered++;
        covered = (HalfSpace**)realloc(covered, numCovered * sizeof(HalfSpace*));
        covered[numCovered - 1] = HalfSpacesMap[hsID];
    } else if (pos == OVERLAPPED) {
        if (isLeaf()) {
            //current node n is a leaf, so we need to insert hs to n's intersectedHS set and then check the capacity of the set
            halfspaces.push_back(hsID);

            //after insertion, max capacity of the set intersectedHalfspace exceeded, we need to split current node
            if (halfspaces.size() > maxCapacity) {
                // Suddividi il nodo
                splitNode();

                // Redistribuisci gli halfspaces tra i figli
                for (auto oldHsID: halfspaces)
                {
                    for (int j = 0; j < numOfSubdivisions; j++) {                     //to each of the 2^d children of root
                        if (this->children[j]->norm)
                            this->children[j]->appendHalfspace(oldHsID);
                    }
                }
                halfspaces.clear();
            }
        } else {
            // Il nodo non è una foglia, passa l'halfspace ai figli
            for (int j = 0; j < numOfSubdivisions; j++) {                     //to each of the 2^d children of root
                if (this->children[j]->norm)
                    this->children[j]->appendHalfspace(hsID);
            }
        }
    }
}

HalfSpace** QNode::getTotalCovered(int& totalCovered) const {
    totalCovered = numCovered;  // Inizia con i covered del nodo corrente
    QNode* ref = parent;

    // Aggiungi i covered degli antenati
    while (ref != nullptr) {
        totalCovered += ref->numCovered;
        ref = ref->parent;
    }

    // Se non ci sono covered, restituisci nullptr
    if (totalCovered == 0) {
        return nullptr;
    }

    // Rialloca memoria per contenere tutti i covered (incluso this)
    auto** totalCoveredArray = (HalfSpace**)malloc(totalCovered * sizeof(HalfSpace*));

    int currentIndex = 0;

    // Copia i covered del nodo corrente (this)
    if (numCovered > 0) {
        memcpy(totalCoveredArray, covered, numCovered * sizeof(HalfSpace*));
        currentIndex += numCovered;
    }

    // Copia i covered degli antenati
    ref = parent;
    while (ref != nullptr) {
        for (int i = 0; i < ref->numCovered; i++) {
            totalCoveredArray[currentIndex++] = ref->covered[i];
        }
        ref = ref->parent;
    }

    return totalCoveredArray;
}


void QNode::splitNode() {
    if (!norm) {
        // Do not split if the node is invalid
        return;
    }

    double*** subDivs = genSubdivisions();

    children = (QNode**)malloc(numOfSubdivisions * sizeof(QNode*));
    // Initialize all children to nullptr
    for (int i = 0; i < numOfSubdivisions; ++i) {
        children[i] = nullptr;
    }

    for (int i = 0; i < numOfSubdivisions; ++i) {
        // Create a new child node
        QNode* child = new QNode(this, subDivs[i], dims);

        // Check if the child node is valid
        if (child->checkNodeValidity()) {
            child->norm = true; // Set the node as valid

        } else {
            // The node is invalid, delete it
            child->norm = false;
        }
        children[i] = child;
        subDivs[i] = nullptr; // The child now owns the mbr
    }

    // Deallocate the subdivisions
    for (int i = 0; i < numOfSubdivisions; ++i) {
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

#define MAXDIM = 10

double*** QNode::genSubdivisions() {
    std::vector<double> mids(dims);
    for (int j = 0; j < dims; ++j) {
        mids[j] = (mbr[j][0] + mbr[j][1]) * 0.5;
    }
    double*** subdivisions = (double***)malloc(numOfSubdivisions * sizeof(double**));

    for (int i = 0; i < numOfSubdivisions; ++i) {
        double** child_mbr = (double**)malloc(dims * sizeof(double*));
        for (int j = 0; j < dims; ++j) {
            child_mbr[j] = (double*)malloc(2 * sizeof(double));
            double mid = mids[j]; // Punto medio dell'MBR nella dimensione j

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
