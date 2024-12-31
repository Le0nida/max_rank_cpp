#include "qtree.h"
#include "qnode.h"
#include "geom.h"
#include "halfspace.h"
#include <iostream>
#include <queue>

extern int numOfSubdivisions;

QTree::QTree(int dims, int maxhsnode)
    : dims(dims), maxhsnode(maxhsnode), root(nullptr)
{
    numOfSubdivisions = (1 << dims); // 2^dims
    root = createroot();
}

QTree::~QTree() {
    // Distruzione dell’intero albero in modo NON ricorsivo
    destroyAllNodes();
    root = nullptr;
}

// Creazione radice
QNode* QTree::createroot() {
    // MBR [0,1]^dims
    std::vector<std::array<double,2>> mbr(dims, {0.0, 1.0});
    QNode* r = new QNode(this, nullptr, mbr);
    return r;
}

void QTree::destroyAllNodes() {
    if (!root) return;
    // Visita BFS (o DFS) di tutti i nodi, e poi delete
    std::queue<QNode*> Q;
    Q.push(root);

    while (!Q.empty()) {
        QNode* current = Q.front();
        Q.pop();

        // Accodiamo i figli (se esistono)
        for (auto child : current->children) {
            if (child) {
                Q.push(child);
            }
        }
        // Eliminiamo il nodo
        delete current;
    }
    // Svuotiamo anche il vettore leaves
    leaves.clear();
}

void QTree::inserthalfspaces(const std::vector<long int>& halfspaces) {
    if (halfspaces.empty() || !root) return;
    root->insertHalfspaces(halfspaces);
}

void QTree::registerLeaf(QNode* leaf) {
    // Inserimento in coda
    leaves.push_back(leaf);
    leaf->leafIndex = (int)leaves.size() - 1;
}

void QTree::unregisterLeaf(QNode* leaf) {
    // Rimozione O(1) tramite swap con l’ultimo e pop_back
    int idx = leaf->leafIndex;
    if (idx >= 0 && idx < (int)leaves.size()) {
        // se leaf non è già l’ultimo, scambiamo
        if (idx != (int)leaves.size()-1) {
            std::swap(leaves[idx], leaves.back());
            // Aggiorniamo l’indice della foglia che era in fondo
            leaves[idx]->leafIndex = idx;
        }
        leaves.pop_back();
        leaf->leafIndex = -1;
    }
}

void QTree::updateAllOrders() {
    if (!root) return;

    root->setOrder(root->covered.size());

    std::queue<QNode*> Q;
    Q.push(root);

    while (!Q.empty()) {
        QNode* current = Q.front();
        Q.pop();

        // figli
        for (auto child : current->children) {
            if (child) {
                size_t newOrder = current->getOrder() + child->covered.size();
                child->setOrder(newOrder);
                Q.push(child);
            }
        }
    }
}
