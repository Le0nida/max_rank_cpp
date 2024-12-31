#ifndef QNODE_H
#define QNODE_H

#include <vector>
#include <array>
#include "halfspace.h"

enum class PositionHS { BELOW, ABOVE, OVERLAPPED };

class QTree;

class QNode {
public:
    QNode(QTree* owner, QNode* parent, const std::vector<std::array<double,2>>& mbr);
    ~QNode();

    QNode(const QNode&) = delete;
    QNode& operator=(const QNode&) = delete;

    bool isRoot() const { return parent == nullptr; }
    bool isLeaf() const { return leaf; }
    void setLeaf(bool lf);

    void setNorm(bool val) { norm = val; }
    bool isNorm() const { return norm; }

    void setOrder(size_t o) { order = o; }
    size_t getOrder() const { return order; }

    const std::vector<std::array<double,2>>& getMBR() const { return mbr; }
    std::vector<long> getHalfspaces() const { return halfspaces; }

    // Ritorna covered locale + (eventuale) risalita del parent se serve
    std::vector<long> getCovered() const;

    // Insert
    void insertHalfspaces(const std::vector<long>& halfspaces);
    void insertHalfspace(long hsID);

    // Splitting
    void splitNode();
    bool checkNodeValidity();

    // MBR vs halfspace
    PositionHS MbrVersusHalfSpace(const std::vector<double>& coeff, double known);

    // Svuota halfspaces
    void clearHalfspaces();

    // Riferimenti
    QTree* owner;
    QNode* parent;

    // Elenco di figli: dimensione 2^dims (alcuni possono essere null)
    // Usiamo un vettore (non dynamic array via new[]).
    // Se preferisci, potresti usare std::array<QNode*, 16> se dims <= 4.
    std::vector<QNode*> children;

    // "covered" e "halfspaces"
    std::vector<long> covered;
    std::vector<long> halfspaces;

    // Indice di questa foglia nel vettore leaves di QTree (se leaf==true)
    int leafIndex;

private:
    std::vector<std::array<double,2>> mbr;
    bool norm;
    bool leaf;
    size_t order;
};

#endif // QNODE_H
