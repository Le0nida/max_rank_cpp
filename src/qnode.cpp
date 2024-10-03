//
// Created by leona on 01/08/2024.
//
#include "qnode.h"
#include "halfspace.h"
#include <iostream>
#include <fstream>
#include <cstring>
#include <numeric>
#include "geom.h"
#include "qtree.h"
#include "qnode.h"
#include "lrucache.h"

int normalizedMax = 1;
int maxCapacity = 10;

QNode::QNode(Context& ctx, int parentID, const std::vector<std::array<double, 2>>& mbr)
    : ctx(ctx), mbr(mbr), leaf(true), norm(true), order(0), parentID(parentID), nodeID(ctx.globalNodeID++) {}

QNode::QNode(Context& ctx)
    : ctx(ctx), mbr({}), leaf(true), norm(true), order(0), parentID(-1), nodeID(ctx.globalNodeID++) {}

QNode::QNode(QNode&& other) noexcept
    : ctx(other.ctx), nodeID(other.nodeID), parentID(other.parentID), mbr(std::move(other.mbr)),
      norm(other.norm), leaf(other.leaf), order(other.order), childrenIDs(std::move(other.childrenIDs)),
      covered(std::move(other.covered)), halfspaces(std::move(other.halfspaces)) {}

QNode& QNode::operator=(QNode&& other) noexcept {
    if (this != &other) {
        nodeID = other.nodeID;
        parentID = other.parentID;
        mbr = std::move(other.mbr);
        norm = other.norm;
        leaf = other.leaf;
        order = other.order;
        childrenIDs = std::move(other.childrenIDs);
        covered = std::move(other.covered);
        halfspaces = std::move(other.halfspaces);
    }
    return *this;
}

int QNode::getNodeID() const {
    return nodeID;
}

int QNode::getParentID() const {
    return parentID;
}

bool QNode::isRoot() const {
    return parentID == -1;
}

bool QNode::isLeaf() const {
    return leaf || childrenIDs.empty();
}

void QNode::setOrder() const {
    size_t localOrder = covered.size();
    std::shared_ptr<QNode> ref = ctx.cache->get(parentID);

    while (ref && ref->parentID >= 0) {
        localOrder += ref->covered.size();
        ref = ctx.cache->get(ref->parentID);
    }

    order = localOrder;
}

std::vector<long int> QNode::getCovered() const {
    std::vector<long int> coveredSpaces = covered;
    std::shared_ptr<QNode> ref = ctx.cache->get(parentID);

    while (ref && !ref->isRoot()) {
        coveredSpaces.insert(coveredSpaces.end(), ref->covered.begin(), ref->covered.end());
        ref = ctx.cache->get(ref->parentID);
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
void QNode::insertHalfspaces(Context& ctx, const std::vector<long int>& halfspaces) {
    size_t num_halfspaces = halfspaces.size();
    for (size_t i = 0; i < num_halfspaces; ++i) {
        long int hsID = halfspaces[i];
        auto hs = ctx.halfspaceCache->get(hsID);
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
                    splitNode(ctx);

                    // Redistribute the intersected halfspaces among the children
                    for (int childID : childrenIDs) {
                        auto child = ctx.cache->get(childID);
                        child->insertHalfspaces(ctx, this->halfspaces);
                    }
                    // Clear the intersected halfspaces from this node
                    this->halfspaces.clear();
                }
            } else {
                // Node is not a leaf, pass the halfspace to children
                for (int childID : childrenIDs) {
                    auto child = ctx.cache->get(childID);
                    child->insertHalfspaces(ctx, {hsID}); // Passing single halfspace as vector
                }
            }
        }/* else if (pos == PositionHS::ABOVE) {
            // The halfspace lies entirely above the MBR, do nothing
            continue;
        }*/
    }
}

// Split the given node
void QNode::splitNode(Context& ctx) {
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
        auto child = std::make_shared<QNode>(ctx);
        child->nodeID = ctx.globalNodeID++;
        //child->level = level + 1;
        child->parentID = nodeID;
        child->mbr = child_mbr;

        // Check if the child node is valid
        if (child->checkNodeValidity()) {
            child->setNorm(true); // Set node as valid
            childrenIDs.push_back(child->getNodeID()); // Add child ID to list
            ctx.cache->add(child); // Add child to cache
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

size_t QNode::saveToDisk(char* pFileData, size_t offset) {
    // Prepare a buffer to hold all data
    std::vector<char> buffer;
    buffer.reserve(estimatedSize()); // Preallocate memory

    // Helper lambda to write data to buffer
    auto writeToBuffer = [&](const void* data, size_t size) {
        const char* bytes = static_cast<const char*>(data);
        buffer.insert(buffer.end(), bytes, bytes + size);
    };

    // Serialize the nodeID and parentID
    writeToBuffer(&nodeID, sizeof(nodeID));
    writeToBuffer(&parentID, sizeof(parentID));

    // Serialize the MBR
    size_t mbrSize = mbr.size();
    writeToBuffer(&mbrSize, sizeof(mbrSize));
    writeToBuffer(mbr.data(), mbrSize * sizeof(std::array<double, 2>));

    // Serialize the norm flag
    writeToBuffer(&norm, sizeof(norm));

    // Serialize the children IDs
    size_t childrenCount = childrenIDs.size();
    writeToBuffer(&childrenCount, sizeof(childrenCount));
    writeToBuffer(childrenIDs.data(), childrenCount * sizeof(long int));

    // Serialize the halfspaces
    size_t halfspaceCount = halfspaces.size();
    writeToBuffer(&halfspaceCount, sizeof(halfspaceCount));
    writeToBuffer(halfspaces.data(), halfspaceCount * sizeof(long int));

    // Serialize the covered halfspaces
    size_t coveredCount = covered.size();
    writeToBuffer(&coveredCount, sizeof(coveredCount));
    writeToBuffer(covered.data(), coveredCount * sizeof(long int));

    // Get total size of the data (excluding the size field)
    size_t dataSize = buffer.size();

    // Write the size of the data
    size_t totalSize = sizeof(dataSize) + dataSize;

    // Copy the data into the memory-mapped file
    char* dest = pFileData + offset;

    // First, write the dataSize
    memcpy(dest, &dataSize, sizeof(dataSize));
    dest += sizeof(dataSize);

    // Then, write the data
    memcpy(dest, buffer.data(), dataSize);

    return totalSize;
}

void QNode::loadFromDisk(char* pFileData, size_t offset) {
    // Read from the memory-mapped file
    const char* src = pFileData + offset;

    // Read the size of the data
    size_t dataSize;
    memcpy(&dataSize, src, sizeof(dataSize));
    src += sizeof(dataSize);

    // Read the data into a buffer
    std::vector<char> buffer(dataSize);
    memcpy(buffer.data(), src, dataSize);

    size_t pos = 0;

    // Helper lambda to read data from buffer
    auto readFromBuffer = [&](void* data, size_t size) {
        memcpy(data, buffer.data() + pos, size);
        pos += size;
    };

    // Deserialize the nodeID and parentID
    readFromBuffer(&nodeID, sizeof(nodeID));
    readFromBuffer(&parentID, sizeof(parentID));

    // Deserialize the MBR
    size_t mbrSize;
    readFromBuffer(&mbrSize, sizeof(mbrSize));
    mbr.resize(mbrSize);
    readFromBuffer(mbr.data(), mbrSize * sizeof(std::array<double, 2>));

    // Deserialize the norm flag
    readFromBuffer(&norm, sizeof(norm));

    // Deserialize the children IDs
    size_t childrenCount;
    readFromBuffer(&childrenCount, sizeof(childrenCount));
    childrenIDs.resize(childrenCount);
    readFromBuffer(childrenIDs.data(), childrenCount * sizeof(long int));

    // Deserialize the halfspaces
    size_t halfspaceCount;
    readFromBuffer(&halfspaceCount, sizeof(halfspaceCount));
    halfspaces.resize(halfspaceCount);
    readFromBuffer(halfspaces.data(), halfspaceCount * sizeof(long int));

    // Deserialize the covered halfspaces
    size_t coveredCount;
    readFromBuffer(&coveredCount, sizeof(coveredCount));
    covered.resize(coveredCount);
    readFromBuffer(covered.data(), coveredCount * sizeof(long int));
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

void QNode::clearHalfspaces() {
    halfspaces.clear();
}

const std::vector<std::array<double, 2>>& QNode::getMBR() const {
    return mbr;
}

bool QNode::isNorm() const {
    return norm;
}

void QNode::setNorm(bool norm) {
    this->norm = norm;
}

void QNode::setLeaf(bool leaf) {
    this->leaf = leaf;
}

size_t QNode::getOrder() const {
    return order;
}

const std::vector<long int>& QNode::getChildrenIDs() const {
    return childrenIDs;
}

void QNode::addChildID(int childID) {
    childrenIDs.push_back(childID);
}

const std::vector<long int>& QNode::getHalfspaces() const {
    return halfspaces;
}

void QNode::setHalfspaces(std::vector<long int> halfspaces) {
    this->halfspaces = std::move(halfspaces);
}
