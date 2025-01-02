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

    // *** NUOVO *** Inserimento parallelo via macro-split:
    // Costruisce una "foresta" di sub-root, uno per ciascuna macro-porzione dell'MBR.
    // Non usa lock: ogni sub-root è costruito in parallelo con i propri halfspaces.
    void inserthalfspacesMacroSplit(const std::vector<long int>& halfspaces);

    // Ritorna TUTTE le foglie (di tutti i sub-root, se stai usando macro-split).
    // Se stai usando la vecchia root, puoi unire i due meccanismi, o ignorare root.
    std::vector<QNode*> getAllLeaves() const;

    // Mantieni se serve
    const std::vector<QNode*>& getLeaves() const { return leaves; }

    // Aggiorna l'ordine di tutti i nodi: sia root che macroRoots
    void updateAllOrders();

    QNode* getRoot() const { return root; }
    int getDims() const { return dims; }
    int getMaxHSNode() const { return maxhsnode; }

private:
    QNode* createroot();
    void destroyAllNodes();

    // *** NUOVO *** Genera i sub-MBR
    std::vector< std::vector<std::array<double,2>> >
    macroSplitMBR(const std::vector<std::array<double,2>>& globalMBR) const;

    // *** NUOVO *** Costruisce in parallelo un sub-root con i halfspaces
    // e restituisce il pointer a QNode*.
    // Usa la logica incrementale interna (insertHalfspaces).
    QNode* buildSubtree(const std::vector<std::array<double,2>>& subMBR,
                        const std::vector<long>& subHS) const;

private:
    int dims;       // dimensioni
    int maxhsnode;  // soglia max halfspaces per splitting

    // Radice old
    QNode* root;

    // Vettore di foglie (approccio "root" old)
    std::vector<QNode*> leaves;

    // *** NUOVO *** "foresta" di sub-root se si usa macroSplit
    std::vector<QNode*> macroRoots;
};

#endif // QTREE_H
