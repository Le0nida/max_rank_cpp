#ifndef QTREE_H
#define QTREE_H

#include "qnode.h"
#include <vector>
#include <queue>

// Numero di partizioni in splitNode (2^dims) dichiarato esternamente
extern int numOfSubdivisions;

class QTree {
public:
    QTree(int dims, int maxhsnode);
    ~QTree();

    // Inserisce un batch di halfspaces
    void inserthalfspaces(const std::vector<long int>& halfspaces);

    // Restituisce (in O(1)) il vettore delle foglie (se vuoi modificare/ordinare: copia o reference)
    const std::vector<QNode*>& getLeaves() const { return leaves; }

    // Registrazione/deregistrazione di foglie (ora su std::vector)
    void registerLeaf(QNode* leaf);
    void unregisterLeaf(QNode* leaf);

    // Aggiorna l'ordine di tutti i nodi
    void updateAllOrders();

    QNode* getRoot() const { return root; }

    int getDims() const { return dims; }
    int getMaxHSNode() const { return maxhsnode; }

private:
    QNode* createroot();

    // Distruzione non-ricorsiva di TUTTI i nodi (BFS/DFS) per evitare stack recursion
    void destroyAllNodes();

    int dims;       // dimensioni
    int maxhsnode;  // soglia max halfspaces per splitting

    // Radice
    QNode* root;

    // Vettore di foglie
    std::vector<QNode*> leaves;
};

#endif // QTREE_H
