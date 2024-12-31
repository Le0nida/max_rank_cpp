#ifndef QNODE_H
#define QNODE_H

#include <vector>
#include <array>
#include <list>
#include "halfspace.h"

enum class PositionHS { BELOW, ABOVE, OVERLAPPED };

// Forward declaration
class QTree;

class QNode {
public:
    // Costruttore: richiede QTree* per potersi registrare come foglia
    QNode(QTree* owner, QNode* parent, const std::vector<std::array<double, 2>>& mbr);
    ~QNode();

    // Disabilita la copia
    QNode(const QNode&) = delete;
    QNode& operator=(const QNode&) = delete;

    // Controlli
    bool isRoot() const { return parent == nullptr; }
    bool isLeaf() const { return leaf; }
    void setLeaf(bool lf);

    // Ordine
    void setOrder(size_t o) { order = o; }
    size_t getOrder() const { return order; }

    // Insert di batch e singolo halfspace
    void insertHalfspaces(const std::vector<long int>& halfspaces);
    void insertHalfspace(long hsID);

    // MBR vs halfspace
    PositionHS MbrVersusHalfSpace(const std::vector<double>& hs_coeff, double hs_known);

    // Splitting del nodo
    void splitNode();

    // Check validity del MBR
    bool checkNodeValidity();

    // Getter vari
    const std::vector<std::array<double, 2>>& getMBR() const { return mbr; }
    bool isNorm() const { return norm; }
    void setNorm(bool val) { norm = val; }

    // Ritorna gli halfspace "coperti" dal nodo (più quelli ereditati dal padre)
    // Usata se serve altrove, ma NON usata in updateAllOrders per evitare overhead
    std::vector<long> getCovered() const;

    std::vector<long> getHalfspaces() const { return halfspaces; }
    // Utility
    void clearHalfspaces();

    // Puntatori e dati pubblici di servizio
    QNode** children;  // array di figli di dimensione numOfSubdivisions
    std::vector<long> covered;      // halfspace completamente "sotto" MBR
    std::vector<long> halfspaces;   // halfspace che effettivamente interessano il nodo

    // Per rimuovere il nodo dalla lista di foglie in O(1)
    std::list<QNode*>::iterator leafIterator;

private:
    // Riferimento all'albero e al parent
    QTree* owner;
    QNode* parent;

    std::vector<std::array<double, 2>> mbr;  // bounding region
    bool norm;   // se il nodo è valido
    bool leaf;   // se è foglia
    size_t order; // "ordine" del nodo

    // Funzioni di supporto
    std::vector<std::vector<std::array<double, 2>>> genSubdivisions();
};

#endif // QNODE_H
