#ifndef QTREE_H
#define QTREE_H

#include "qnode.h"
#include "halfspace.h"
#include <vector>
#include <array>
#include <list>
#include <queue>

// Numero di partizioni in splitNode() = 2^dims (es. in un octree 2^3=8).
extern int numOfSubdivisions;

class QTree {
public:
    QTree(int dims, int maxhsnode);
    ~QTree();

    // Inserisce un batch di halfspaces nella radice (e ricorsivamente nei figli)
    void inserthalfspaces(const std::vector<long int>& halfspaces);

    // Restituisce la lista delle foglie aggiornata in O(1)
    std::list<QNode*> getLeaves() const { return leaves; }

    // Metodi per gestire la registrazione/deregistrazione di un nodo-foglia
    void registerLeaf(QNode* leaf);
    void unregisterLeaf(QNode* leaf);

    // Funzione per aggiornare l'ordine di tutti i nodi in BFS
    void updateAllOrders();

    // Getter sul root
    QNode* getRoot() const { return root; }

    // Parametri
    int getDims() const { return dims; }
    int getMaxHSNode() const { return maxhsnode; }

private:
    QNode* createroot();

    int dims;       // Dimensioni dello spazio
    int maxhsnode;  // Capacità massima halfspaces in un nodo

    // Radice dell'albero
    QNode* root;

    // Lista aggiornata in tempo reale di foglie
    std::list<QNode*> leaves;
};

#endif // QTREE_H
