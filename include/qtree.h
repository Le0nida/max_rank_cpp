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

    // Inserisce in modo tradizionale i halfspaces (incrementale)
    void inserthalfspaces(const std::vector<long int>& halfspaces);

    // Inserimento parallelo via macro-split
    void inserthalfspacesMacroSplit(const std::vector<long int>& halfspaces);

    // Ritorna TUTTE le foglie (di tutti i sub-root)
    std::vector<QNode*> getAllLeaves() const;

    // Ritorna le foglie della "old root" (se usata)
    const std::vector<QNode*>& getLeaves() const { return leaves; }

    // Aggiorna l'ordine di tutti i nodi: root + macroRoots
    void updateAllOrders();

    QNode* getRoot() const { return root; }
    int getDims() const { return dims; }
    int getMaxHSNode() const { return maxhsnode; }

private:
    QNode* createroot();
    void destroyAllNodes();

    // Genera (una volta sola) i sub-MBR e li memorizza in precomputedSubMBRs
    std::vector< std::vector<std::array<double,2>> >
    macroSplitMBR(const std::vector<std::array<double,2>>& globalMBR) const;

    // Costruisce un sub-root con i halfspaces
    QNode* buildSubtree(const std::vector<std::array<double,2>>& subMBR,
                        const std::vector<long>& subHS) const;

private:
    int dims;       // dimensioni
    int maxhsnode;  // soglia max halfspaces per splitting

    // Radice "classica" (eventualmente usata)
    QNode* root;

    // Vettore di foglie relativo alla root classica
    std::vector<QNode*> leaves;

    // "foresta" di sub-root se si usa macroSplit
    std::vector<QNode*> macroRoots;

    // **NOVITÀ**: SubMBR pre-calcolati (dimensione 2^dims)
    std::vector< std::vector<std::array<double,2>> > precomputedSubMBRs;
};

#endif // QTREE_H
