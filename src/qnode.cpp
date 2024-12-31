#include "qnode.h"
#include "qtree.h"
#include "geom.h"
#include <algorithm>
#include <iostream>
#include <queue>

// Valore di riferimento per la checkNodeValidity()
static int normalizedMax = 1;
// Capacità massima del nodo, al superamento si fa split
static int maxCapacity = 10;

extern int numOfSubdivisions;

// Contatore globale per ID (se necessario)
static int globalNodeID = 0;

// Costruttore: registra immediatamente il nodo come foglia
QNode::QNode(QTree* owner, QNode* parent, const std::vector<std::array<double, 2>>& mbr)
  : owner(owner),
    parent(parent),
    mbr(mbr),
    norm(true),
    leaf(true),
    order(0)
{
    // Alloca i figli
    children = nullptr;
    if (owner) {
        // Ci registriamo come foglia presso l'albero
        owner->registerLeaf(this);
    }
    // Se la dimensione halfspaces del nodo dovesse superare maxCapacity,
    // eseguiamo splitNode() dinamicamente in insertHalfspace().
}

// Distruttore ricorsivo
QNode::~QNode() {
    // Se è foglia, deregistra dalla lista delle foglie
    if (leaf && owner) {
        owner->unregisterLeaf(this);
    }
    // Dealloca i figli
    if (children) {
        for (int i = 0; i < numOfSubdivisions; ++i) {
            if (children[i]) {
                delete children[i];
            }
        }
        delete[] children;
        children = nullptr;
    }
}

// setLeaf con registrazione/deregistrazione come nella “prima versione”
void QNode::setLeaf(bool lf) {
    if (leaf == lf) return; // Nessun cambiamento
    leaf = lf;
    if (leaf) {
        // registrazione nella lista di foglie
        if (owner) owner->registerLeaf(this);
    } else {
        // deregistrazione dalla lista di foglie
        if (owner) owner->unregisterLeaf(this);
    }
}

PositionHS QNode::MbrVersusHalfSpace(const std::vector<double>& hs_coeff, double hs_known) {
    if (mbr.empty()) {
        return PositionHS::OVERLAPPED;
    }

    double minVal = 0.0;
    double maxVal = 0.0;
    size_t dims = hs_coeff.size();

    for (size_t i = 0; i < dims; ++i) {
        double coeff = hs_coeff[i];
        if (coeff >= 0) {
            minVal += coeff * mbr[i][0];
            maxVal += coeff * mbr[i][1];
        } else {
            minVal += coeff * mbr[i][1];
            maxVal += coeff * mbr[i][0];
        }
    }

    if (maxVal < hs_known) {
        return PositionHS::BELOW;
    } else if (minVal > hs_known) {
        return PositionHS::ABOVE;
    } else {
        return PositionHS::OVERLAPPED;
    }
}

// Inserimento di un singolo halfspace
void QNode::insertHalfspace(long hsID) {
    auto hs = halfspaceCache->get(hsID);
    if (!hs) return; // Non esiste

    PositionHS pos = MbrVersusHalfSpace(hs->coeff, hs->known);

    switch (pos) {
    case PositionHS::BELOW:
        covered.push_back(hsID);
        break;
    case PositionHS::OVERLAPPED:
        if (isLeaf()) {
            // Aggiungiamo agli halfspaces del nodo
            halfspaces.push_back(hsID);

            // Se sforiamo la capacità, splittiamo
            if ((int)halfspaces.size() > maxCapacity) {
                splitNode();
                // Dopo lo split, se i figli sono validi, ridistribuiamo
                if (children) {
                    for (auto h : halfspaces) {
                        for (int k = 0; k < numOfSubdivisions; k++) {
                            if (children[k]) {
                                children[k]->insertHalfspace(h);
                            }
                        }
                    }
                    halfspaces.clear();
                    halfspaces.shrink_to_fit();
                }
            }
        } else {
            // Propaghiamo ai figli
            if (children) {
                for (int k = 0; k < numOfSubdivisions; k++) {
                    if (children[k]) {
                        children[k]->insertHalfspace(hsID);
                    }
                }
            }
        }
        break;
    case PositionHS::ABOVE:
        // Non interessa questo MBR
        break;
    }
}

// Inserimento di un batch di halfspaces (evitiamo vector temporanei)
void QNode::insertHalfspaces(const std::vector<long int>& new_halfspaces) {
    for (auto hsID : new_halfspaces) {
        insertHalfspace(hsID);
    }
}

void QNode::splitNode() {
    if (!norm) return;

    // Se abbiamo pochi halfspace, non splitto
    size_t totalHS = halfspaces.size() + covered.size();
    if (totalHS < (size_t)maxCapacity) {
        return;
    }

    // Non siamo più foglia
    setLeaf(false);

    // Generiamo i MBR dei figli
    std::vector<std::vector<std::array<double, 2>>> subDivs = genSubdivisions();
    children = new QNode*[numOfSubdivisions];
    for (int i = 0; i < numOfSubdivisions; ++i) {
        children[i] = nullptr;
    }

    int i = 0;
    for (auto& child_mbr : subDivs) {
        QNode* child = new QNode(owner, this, child_mbr);
        if (child->checkNodeValidity()) {
            child->setNorm(true);
            children[i] = child;
        } else {
            // nodo invalido, lo distruggiamo
            child->setNorm(false);
            delete child;
            children[i] = nullptr;
        }
        i++;
    }
}

// Genera 2^dims subdiv per splittare l'MBR
std::vector<std::vector<std::array<double,2>>> QNode::genSubdivisions() {
    size_t dims = mbr.size();
    std::vector<std::vector<std::array<double,2>>> subdivisions;
    subdivisions.reserve(numOfSubdivisions); // evitiamo un po' di realloc

    for (int mask = 0; mask < numOfSubdivisions; ++mask) {
        std::vector<std::array<double,2>> child_mbr(dims);
        for (size_t d = 0; d < dims; d++) {
            double mid = (mbr[d][0] + mbr[d][1]) * 0.5;
            if (mask & (1 << d)) {
                // metà superiore
                child_mbr[d] = {mid, mbr[d][1]};
            } else {
                // metà inferiore
                child_mbr[d] = {mbr[d][0], mid};
            }
        }
        subdivisions.push_back(child_mbr);
    }
    return subdivisions;
}

// Controlla se almeno un vertice dell'MBR soddisfa sum <= normalizedMax
bool QNode::checkNodeValidity() {
    size_t dims = mbr.size();
    size_t num_vertices = (1 << dims);

    for (size_t i = 0; i < num_vertices; ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < dims; j++) {
            double coord = (i & (1 << j)) ? mbr[j][1] : mbr[j][0];
            sum += coord;
        }
        if (sum <= normalizedMax) {
            return true; // almeno un vertice valido
        }
    }
    return false;
}

// Restituisce covered "locale" + quelli dei genitori (se si desidera).
// ATTENZIONE: NON viene usato in updateAllOrders, per evitare eccessive copie
std::vector<long> QNode::getCovered() const {
    std::vector<long> out(covered.begin(), covered.end());

    // Se vogliamo aggiungere i covered di tutti i parent
    const QNode* cur = parent;
    while (cur) {
        out.insert(out.end(), cur->covered.begin(), cur->covered.end());
        cur = cur->parent;
    }
    return out;
}

void QNode::clearHalfspaces() {
    halfspaces.clear();
    halfspaces.shrink_to_fit();
}
