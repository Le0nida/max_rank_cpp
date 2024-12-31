#include "qnode.h"
#include "qtree.h"
#include "geom.h"
#include <iostream>
#include <algorithm>
#include <queue>

static int normalizedMax = 1;
static int maxCapacity = 10;  // <= potresti volerlo più alto!
extern int numOfSubdivisions;

QNode::QNode(QTree* owner, QNode* parent, const std::vector<std::array<double,2>>& mbr)
    : owner(owner),
      parent(parent),
      children(), // vuoto, verrà riempito in splitNode() se serve
      leafIndex(-1),
      mbr(mbr),
      norm(true),
      leaf(true),
      order(0)
{
    if (owner) {
        owner->registerLeaf(this); // Mi registro come foglia
    }
}

QNode::~QNode() {
    // ATTENZIONE: Non eliminiamo i figli qui, perché lo fa QTree in modo non ricorsivo.
    // Se siamo foglia, deregistriamo
    if (leaf && owner) {
        owner->unregisterLeaf(this);
    }
}

void QNode::setLeaf(bool lf) {
    if (leaf == lf) return;
    leaf = lf;
    if (leaf) {
        if (owner) owner->registerLeaf(this);
    } else {
        if (owner) owner->unregisterLeaf(this);
    }
}

std::vector<long> QNode::getCovered() const {
    // Se serve aggregare i covered dei parent:
    std::vector<long> out(covered.begin(), covered.end());
    const QNode* anc = parent;
    while (anc) {
        out.insert(out.end(), anc->covered.begin(), anc->covered.end());
        anc = anc->parent;
    }
    return out;
}

PositionHS QNode::MbrVersusHalfSpace(const std::vector<double>& coeff, double known) {
    if (mbr.empty()) {
        return PositionHS::OVERLAPPED;
    }
    double minVal = 0.0;
    double maxVal = 0.0;
    size_t dims = coeff.size();

    for (size_t i = 0; i < dims; ++i) {
        double c = coeff[i];
        if (c >= 0) {
            minVal += c * mbr[i][0];
            maxVal += c * mbr[i][1];
        } else {
            minVal += c * mbr[i][1];
            maxVal += c * mbr[i][0];
        }
    }
    if (maxVal < known) return PositionHS::BELOW;
    if (minVal > known) return PositionHS::ABOVE;
    return PositionHS::OVERLAPPED;
}

void QNode::insertHalfspace(long hsID) {
    auto hs = halfspaceCache->get(hsID);
    if (!hs) return;

    PositionHS pos = MbrVersusHalfSpace(hs->coeff, hs->known);

    switch (pos) {
        case PositionHS::BELOW:
            covered.push_back(hsID);
            break;
        case PositionHS::OVERLAPPED:
            if (isLeaf()) {
                halfspaces.push_back(hsID);
                // Se sforiamo la capacità => split
                if ((int)halfspaces.size() > maxCapacity && norm) {
                    splitNode();
                    // Se abbiamo figli, ridistribuiamo
                    if (!children.empty()) {
                        for (auto h : halfspaces) {
                            for (auto* ch : children) {
                                if (ch) ch->insertHalfspace(h);
                            }
                        }
                        halfspaces.clear();
                        halfspaces.shrink_to_fit();
                    }
                }
            } else {
                // Propaghiamo
                for (auto* ch : children) {
                    if (ch) ch->insertHalfspace(hsID);
                }
            }
            break;
        case PositionHS::ABOVE:
            // Non lo gestiamo
            break;
    }
}

void QNode::insertHalfspaces(const std::vector<long>& new_halfspaces) {
    for (auto id : new_halfspaces) {
        insertHalfspace(id);
    }
}

void QNode::splitNode() {
    // Già controllato da caller se n>maxCapacity, ma ricontrolliamo
    if (!norm) return;

    size_t totalHS = halfspaces.size() + covered.size();
    if (totalHS < (size_t)maxCapacity) {
        return;
    }

    // Diventiamo non-foglia
    setLeaf(false);

    // Preallochiamo vector di figli
    children.resize(numOfSubdivisions, nullptr);

    // Generiamo i MBR
    size_t dims = mbr.size();
    for (int mask = 0; mask < numOfSubdivisions; ++mask) {
        std::vector<std::array<double,2>> child_mbr(dims);
        for (size_t d = 0; d < dims; d++) {
            double mid = 0.5*(mbr[d][0] + mbr[d][1]);
            if (mask & (1 << d)) {
                child_mbr[d] = {mid, mbr[d][1]};
            } else {
                child_mbr[d] = {mbr[d][0], mid};
            }
        }
        // Creo QNode figlio
        QNode* child = new QNode(owner, this, child_mbr);
        if (child->checkNodeValidity()) {
            child->setNorm(true);
            children[mask] = child;
        } else {
            child->setNorm(false);
            delete child; // niente ricorsione, è un singolo new
            children[mask] = nullptr;
        }
    }
}

bool QNode::checkNodeValidity() {
    size_t dims = mbr.size();
    size_t num_vertices = (1 << dims);

    for (size_t i = 0; i < num_vertices; ++i) {
        double sum = 0.0;
        for (size_t d = 0; d < dims; d++) {
            double coord = (i & (1 << d)) ? mbr[d][1] : mbr[d][0];
            sum += coord;
        }
        if (sum <= normalizedMax) {
            return true;
        }
    }
    return false;
}

void QNode::clearHalfspaces() {
    halfspaces.clear();
    halfspaces.shrink_to_fit();
}
