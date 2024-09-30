//
// Created by leona on 01/08/2024.
//

#include "qnode.h"
#include "halfspace.h"
#include <iostream>
#include <fstream>
#include <cstring>  // Per memcpy
#include "lrucache.h"
#include <numeric>
#include "geom.h"
#include "qtree.h"

int normalizedMax = 1;
int maxCapacity = 10; // sarebbe maxhsnode

// Counter globale per assegnare ID unici (opzione semplice)
int globalNodeID = 0;

// Constructor for QNode, initializes the member variables.
QNode::QNode(int parentID, const std::vector<std::array<double, 2>>& mbr)
    : mbr(mbr), leaf(true), norm(true), order(0), parentID(parentID), nodeID(globalNodeID++){}

QNode::QNode()
    : mbr({}), leaf(true), norm(true), order(0), parentID(0), nodeID(-1){}


// Checks if the node is the root of the tree.
bool QNode::isRoot() const {
    return parentID == -1;
}

// Checks if the node is a leaf (i.e., has no children).
bool QNode::isLeaf() const {
    return leaf || childrenIDs.empty();
}

// Computes the order of the node by traversing back up the tree.
void QNode::setOrder() const {
    size_t localOrder = covered.size();  // Initialize order with the size of covered halfspaces in the current node.
    std::shared_ptr<QNode> ref = globalCache.get(parentID);

    // Traverse back up the tree to accumulate the order from parent nodes.
    while (ref->parentID >= 0) {
        localOrder += ref->covered.size();
        ref = globalCache.get(ref->parentID);
    }

    order = localOrder;
}

void QNode::clearHalfspaces()
{
    if (!halfspaces.empty()) {
        halfspaces.clear();
    }
}


// Retrieves the covering halfspaces by traversing back up the tree.
std::vector<long int> QNode::getCovered() const{
    std::vector<long int> coveredSpaces = covered;  // Start with halfspaces covered in the current node.
    std::shared_ptr<QNode> ref = globalCache.get(parentID);;

    // Traverse back up the tree to accumulate the covered halfspaces from parent nodes.
    while (ref && !ref->isRoot()) {
        coveredSpaces.insert(coveredSpaces.end(), ref->covered.begin(), ref->covered.end());
        ref = globalCache.get(ref->parentID);;
    }

    return coveredSpaces;
}

PositionHS QNode::MbrVersusHalfSpace(const std::vector<double>& hs_coeff, double hs_known, const std::vector<std::array<double, 2>>& mbr) {
    double minVal = 0.0;
    double maxVal = 0.0;
    size_t dims = hs_coeff.size();

    // Compute the minimum and maximum values of the halfspace over the MBR
    for (size_t i = 0; i < dims; ++i) {
        double coeff = hs_coeff[i];
        if (coeff >= 0) {
            minVal += coeff * mbr[i][0]; // Use minimum MBR value for positive coefficients
            maxVal += coeff * mbr[i][1]; // Use maximum MBR value for positive coefficients
        } else {
            minVal += coeff * mbr[i][1]; // Use maximum MBR value for negative coefficients
            maxVal += coeff * mbr[i][0]; // Use minimum MBR value for negative coefficients
        }
    }

    // Compare min and max values with hs_known to determine position
    if (maxVal < hs_known) {
        // The entire MBR lies below the halfspace
        return PositionHS::BELOW;
    } else if (minVal > hs_known) {
        // The entire MBR lies above the halfspace
        return PositionHS::ABOVE;
    } else {
        // The halfspace overlaps with the MBR
        return PositionHS::OVERLAPPED;
    }
}

// Inserts halfspaces into the node and distributes them to children nodes if necessary.
void QNode::insertHalfspaces(const std::vector<long int>& halfspaces) {
    size_t num_halfspaces = halfspaces.size();
    for (size_t i = 0; i < num_halfspaces; ++i) {
        long int hsID = halfspaces[i];
        auto hs = halfspaceCache->get(hsID);
        const std::vector<double>& hs_coeff = hs->coeff;
        double hs_known = hs->known;

        // Determine the PositionHS of the halfspace relative to the node's MBR
        PositionHS pos = MbrVersusHalfSpace(hs_coeff, hs_known, mbr);

        if (pos == PositionHS::BELOW) {
            // Halfspace lies entirely below the MBR, store in covered halfspaces
            covered.push_back(hsID);
        } else if (pos == PositionHS::OVERLAPPED) {
            if (isLeaf()) {
                // Node is a leaf, store in intersectedHalfspaces
                this->halfspaces.push_back(hsID);

                // Check if we need to split the node
                if (this->halfspaces.size() > maxCapacity) {
                    // Split the node
                    splitNode();

                    // Redistribute the intersected halfspaces among the children
                    for (int childID : childrenIDs) {
                        auto child = globalCache.get(childID);
                        child->insertHalfspaces(this->halfspaces);
                    }
                    // Clear the intersected halfspaces from this node
                    this->halfspaces.clear();
                }
            } else {
                // Node is not a leaf, pass the halfspace to children
                for (int childID : childrenIDs) {
                    auto child = globalCache.get(childID);
                    child->insertHalfspaces({hsID}); // Passing single halfspace as vector
                }
            }
        }/* else if (pos == PositionHS::ABOVE) {
            // The halfspace lies entirely above the MBR, do nothing
            continue;
        }*/
    }
}

// Split the given node
void QNode::splitNode() {
    /*if (level >= maxLevel) {
        // Do not split if maximum level reached
        return;
    }*/

    if (!isNorm()) {
        // Do not split if the node is invalid
        return;
    }

    // Generate subdivisions (child MBRs)
    std::vector<std::vector<std::array<double, 2>>> subDivs = genSubdivisions();

    childrenIDs.clear(); // Clear any existing children IDs

    for (const auto& child_mbr : subDivs) {
        // Create new child node
        auto child = std::make_shared<QNode>(/* appropriate parameters */);
        child->nodeID = globalNodeID++;
        //child->level = level + 1;
        child->parentID = nodeID;
        child->mbr = child_mbr;

        // Check if the child node is valid
        if (child->checkNodeValidity()) {
            child->setNorm(true); // Set node as valid
            childrenIDs.push_back(child->getNodeID()); // Add child ID to list
            globalCache.add(child); // Add child to cache
        } else {
            // Node is invalid, do not add to children
            continue;
        }
    }

    setLeaf(false); // This node is no longer a leaf
}

std::vector<std::vector<std::array<double, 2>>> QNode::genSubdivisions() {
    size_t dims = mbr.size();
    std::vector<std::vector<std::array<double, 2>>> subdivisions;

    size_t numOfSubdivisions = 1 << dims; // 2^dims combinations

    for (int i = 0; i < numOfSubdivisions; ++i) {
        std::vector<std::array<double, 2>> child_mbr(dims);

        for (int j = 0; j < dims; ++j) {
            double mid = (mbr[j][0] + mbr[j][1]) / 2.0; // Midpoint of MBR in dimension j

            if (i & (1 << j)) {
                // Use upper half in dimension j
                child_mbr[j][0] = mid;
                child_mbr[j][1] = mbr[j][1];
            } else {
                // Use lower half in dimension j
                child_mbr[j][0] = mbr[j][0];
                child_mbr[j][1] = mid;
            }
        }

        subdivisions.push_back(child_mbr);
    }

    return subdivisions;
}

bool QNode::checkNodeValidity() {
    size_t dims = mbr.size();
    size_t num_vertices = 1 << dims; // Number of vertices is 2^dims

    for (size_t i = 0; i < num_vertices; ++i) {
        double sum = 0.0;

        // Generate the coordinates of the vertex
        for (size_t j = 0; j < dims; ++j) {
            double coord = (i & (1 << j)) ? mbr[j][1] : mbr[j][0]; // Use upper or lower bound based on bitmask
            sum += coord;
        }

        if (sum <= normalizedMax) {
            // At least one vertex is valid
            return true;
        }
    }
    // All vertices are invalid
    return false;
}

// Serialize the QNode to the provided output stream (in binary)
void QNode::serialize(std::ostream& outStream) const {
    // Serialize basic fields
    outStream.write(reinterpret_cast<const char*>(&nodeID), sizeof(nodeID));
    outStream.write(reinterpret_cast<const char*>(&parentID), sizeof(parentID));
    outStream.write(reinterpret_cast<const char*>(&norm), sizeof(norm));
    outStream.write(reinterpret_cast<const char*>(&leaf), sizeof(leaf));
    outStream.write(reinterpret_cast<const char*>(&order), sizeof(order));

    // Serialize MBR (Minimum Bounding Region)
    size_t mbrSize = mbr.size();
    outStream.write(reinterpret_cast<const char*>(&mbrSize), sizeof(mbrSize));
    for (const auto& arr : mbr) {
        outStream.write(reinterpret_cast<const char*>(&arr[0]), sizeof(double));
        outStream.write(reinterpret_cast<const char*>(&arr[1]), sizeof(double));
    }

    // Serialize children IDs
    size_t childrenCount = childrenIDs.size();
    outStream.write(reinterpret_cast<const char*>(&childrenCount), sizeof(childrenCount));
    outStream.write(reinterpret_cast<const char*>(childrenIDs.data()), childrenCount * sizeof(long int));

    // Serialize halfspaces
    size_t halfspaceCount = halfspaces.size();
    outStream.write(reinterpret_cast<const char*>(&halfspaceCount), sizeof(halfspaceCount));
    outStream.write(reinterpret_cast<const char*>(halfspaces.data()), halfspaceCount * sizeof(long int));

    // Serialize covered halfspaces
    size_t coveredCount = covered.size();
    outStream.write(reinterpret_cast<const char*>(&coveredCount), sizeof(coveredCount));
    outStream.write(reinterpret_cast<const char*>(covered.data()), coveredCount * sizeof(long int));
}

// Deserialize the QNode from the provided input stream (in binary)
void QNode::deserialize(std::istream& inStream) {
    // Deserialize basic fields
    inStream.read(reinterpret_cast<char*>(&nodeID), sizeof(nodeID));
    inStream.read(reinterpret_cast<char*>(&parentID), sizeof(parentID));
    inStream.read(reinterpret_cast<char*>(&norm), sizeof(norm));
    inStream.read(reinterpret_cast<char*>(&leaf), sizeof(leaf));
    inStream.read(reinterpret_cast<char*>(&order), sizeof(order));

    // Deserialize MBR (Minimum Bounding Region)
    size_t mbrSize;
    inStream.read(reinterpret_cast<char*>(&mbrSize), sizeof(mbrSize));
    mbr.resize(mbrSize);
    for (auto& arr : mbr) {
        inStream.read(reinterpret_cast<char*>(&arr[0]), sizeof(double));
        inStream.read(reinterpret_cast<char*>(&arr[1]), sizeof(double));
    }

    // Deserialize children IDs
    size_t childrenCount;
    inStream.read(reinterpret_cast<char*>(&childrenCount), sizeof(childrenCount));
    childrenIDs.resize(childrenCount);
    inStream.read(reinterpret_cast<char*>(childrenIDs.data()), childrenCount * sizeof(long int));

    // Deserialize halfspaces
    size_t halfspaceCount;
    inStream.read(reinterpret_cast<char*>(&halfspaceCount), sizeof(halfspaceCount));
    halfspaces.resize(halfspaceCount);
    inStream.read(reinterpret_cast<char*>(halfspaces.data()), halfspaceCount * sizeof(long int));

    // Deserialize covered halfspaces
    size_t coveredCount;
    inStream.read(reinterpret_cast<char*>(&coveredCount), sizeof(coveredCount));
    covered.resize(coveredCount);
    inStream.read(reinterpret_cast<char*>(covered.data()), coveredCount * sizeof(long int));
}


size_t QNode::estimatedSize() const {
    size_t size = 0;
    size += sizeof(nodeID);
    size += sizeof(parentID);
    size += sizeof(size_t); // mbrSize
    size += mbr.size() * sizeof(std::array<double, 2>);
    size += sizeof(norm);
    size += sizeof(size_t); // childrenCount
    size += childrenIDs.size() * sizeof(long int);
    size += sizeof(size_t); // halfspaceCount
    size += halfspaces.size() * sizeof(long int);
    size += sizeof(size_t); // coveredCount
    size += covered.size() * sizeof(long int);
    return size;
}