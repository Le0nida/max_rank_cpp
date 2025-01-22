#include "qnode.h"
#include "qtree.h"

static int normalizedMax = 1;

QNode::QNode(QTree* owner,
             QNode* parent,
             const std::vector<std::array<float,2>>& mbr,
             const int level)
    : owner(owner),
      parent(parent),
      children(),
      covered(),
      halfspaces(),
      leafIndex(-1),
      mbr(mbr),
      norm(true),
      leaf(true),
      order(0),
      level(level)
{
}

void QNode::setLeaf(const bool lf) {
    leaf = lf;
}

std::vector<long> QNode::getCovered() const {
    // Combine local covered with ancestor's covered (delta coverage approach)
    std::vector<long> out;
    out.reserve(covered.size() + 32); // slight guess, to reduce reallocation
    out.insert(out.end(), covered.begin(), covered.end());

    const QNode* anc = parent;
    while (anc) {
        out.insert(out.end(), anc->covered.begin(), anc->covered.end());
        anc = anc->parent;
    }
    return out;
}

PositionHS QNode::MbrVersusHalfSpace(const std::vector<double>& coeff, const double known) const {
    if (mbr.empty()) {
        return PositionHS::OVERLAPPED;
    }
    double minVal = 0.0, maxVal = 0.0;
    size_t dcount = coeff.size();

    for (size_t i = 0; i < dcount; ++i) {
        double c = coeff[i];
        // MBR coords in float
        float low  = mbr[i][0];
        float high = mbr[i][1];
        if (c >= 0) {
            minVal += c * low;
            maxVal += c * high;
        } else {
            minVal += c * high;
            maxVal += c * low;
        }
    }
    if (maxVal < known)  return PositionHS::BELOW;
    if (minVal > known)  return PositionHS::ABOVE;
    return PositionHS::OVERLAPPED;
}

void QNode::insertHalfspace(const long hsID) {
    const auto hs = halfspaceCache->get(hsID);
    if (!hs) return;

    PositionHS pos = MbrVersusHalfSpace(hs->coeff, hs->known);
    switch (pos) {
        case PositionHS::BELOW:
            // This halfspace is fully covering the node (delta coverage)
            covered.push_back(hsID);
            break;
        case PositionHS::OVERLAPPED:
            // Partially covers the node
            if (leaf) {
                halfspaces.push_back(hsID);
                // Check capacity -> split if we exceed maxhsnode
                if ((int)halfspaces.size() > owner->maxhsnode && norm) {
                    if (level < owner->maxLevel) {
                        splitNode();
                        // Redistribute existing halfspaces to children
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
                }
            } else {
                // Propagate to children
                for (auto* ch : children) {
                    if (ch) ch->insertHalfspace(hsID);
                }
            }
            break;
        case PositionHS::ABOVE:
            // Node is completely outside the halfspace -> do nothing
            break;
    }
}

void QNode::insertHalfspaces(const std::vector<long>& new_halfspaces) {
    for (const auto id : new_halfspaces) {
        insertHalfspace(id);
    }
}

void QNode::splitNode() {
    // Do not split if at max level or not valid
    if (level == owner->maxLevel || !norm) {
        return;
    }

    // total halfspaces
    size_t totalHS = halfspaces.size() + covered.size();
    if (totalHS < (size_t)owner->maxhsnode) return;

    // Become an internal node
    setLeaf(false);
    children.resize(numOfSubdivisions, nullptr);

    size_t dcount = mbr.size();
    for (int mask = 0; mask < numOfSubdivisions; ++mask) {
        std::vector<std::array<float,2>> child_mbr(dcount);
        for (size_t d = 0; d < dcount; d++) {
            float minVal = mbr[d][0];
            float maxVal = mbr[d][1];
            float mid = 0.5f * (minVal + maxVal);
            if (mask & (1 << d)) {
                child_mbr[d] = { mid, maxVal };
            } else {
                child_mbr[d] = { minVal, mid };
            }
        }
        // Create the child node
        auto* child = new QNode(owner, this, child_mbr, level + 1);

        // Validate the child node
        if (child->checkNodeValidity()) {
            child->norm = true;
            children[mask] = child;
        } else {
            child->norm = false;
            delete child;
            children[mask] = nullptr;
        }
    }
}

bool QNode::checkNodeValidity() const {
    // For each corner in [min,max], check if sum of coordinates <= normalizedMax
    size_t dcount = mbr.size();
    size_t corners = (1 << dcount);

    for (size_t i = 0; i < corners; ++i) {
        double sum = 0.0;
        for (size_t d = 0; d < dcount; d++) {
            float coord = (i & (1 << d)) ? mbr[d][1] : mbr[d][0];
            sum += coord;
        }
        // If at least one corner is inside normalized range, consider it valid
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
