#include "qtree.h"
#include "qnode.h"
#include "geom.h"
#include "halfspace.h"
#include <numeric>
#include <iostream>

// Dichiarazione esterna: quante partizioni si generano in splitNode()
// (es. se dims=2, numOfSubdivisions=4; se dims=3, numOfSubdivisions=8; etc.)


QTree::QTree(int dims, int maxhsnode)
    : dims(dims), maxhsnode(maxhsnode), root(nullptr)
{
    // Calcoliamo 2^dims
    numOfSubdivisions = (1 << dims);

    // Creiamo la radice
    root = createroot();
}

QTree::~QTree() {
    // Distruggiamo tutta la gerarchia, cominciando dalla radice
    // (il distruttore di QNode farà il delete ricorsivo dei figli).
    if (root) delete root;
    root = nullptr;
}

// Creazione della radice e suo eventuale primo split
QNode* QTree::createroot() {
    // MBR [0,1]^dims
    std::vector<std::array<double, 2>> mbr(dims, {0.0, 1.0});

    // Crea QNode radice (parent=nullptr, isLeaf=true di default)
    QNode* r = new QNode(this, nullptr, mbr);
    return r;
}

// Inserisce i nuovi halfspaces (bypassando i costruttori di vector temporanei)
void QTree::inserthalfspaces(const std::vector<long int>& halfspaces) {
    if (halfspaces.empty()) return;
    // Radice fa la insert
    root->insertHalfspaces(halfspaces);
}

// Registra un nuovo nodo foglia
void QTree::registerLeaf(QNode* leaf) {
    // Aggiunge in coda alla lista
    leaf->leafIterator = leaves.insert(leaves.end(), leaf);
}

// Deregistra un nodo che non è più foglia
void QTree::unregisterLeaf(QNode* leaf) {
    // Rimozione O(1) tramite l'iteratore salvato
    if (leaf->leafIterator != leaves.end()) {
        leaves.erase(leaf->leafIterator);
        leaf->leafIterator = leaves.end();
    }
}

// Aggiorna l'ordine di tutti i nodi in BFS, evitando copie di vettori
void QTree::updateAllOrders() {
    if (!root) return;

    // Inizializziamo l'ordine del root come la dimensione di covered
    root->setOrder(root->covered.size());

    // BFS con coda
    std::queue<QNode*> Q;
    Q.push(root);

    while (!Q.empty()) {
        QNode* current = Q.front();
        Q.pop();

        // Per ogni figlio
        if (current->children)
        {
            for (int i = 0; i < numOfSubdivisions; ++i) {
                QNode* child = current->children[i];
                if (!child) continue;

                // L'ordine del figlio = ordine del padre + # covered locali del figlio
                child->setOrder(current->getOrder() + child->covered.size());
                Q.push(child);
            }
        }
    }
}
