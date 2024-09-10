//
// Created by leona on 01/08/2024.
//

#ifndef QNODE_H
#define QNODE_H

#include <vector>
#include <array>
#include <memory>
#include "halfspace.h"

class QNode {
public:
    QNode(int parentID, const std::vector<std::array<double, 2>>& mbr);
    ~QNode() = default;
    explicit QNode(int id) : nodeID(id) {}
    QNode() = default;

    QNode(const QNode&) = delete;               // Disabilita il costruttore di copia
    QNode& operator=(const QNode&) = delete;    // Disabilita l'assegnazione tramite copia


    QNode(QNode&& other) noexcept
    : nodeID(other.nodeID),
      parentID(other.parentID),
      mbr(std::move(other.mbr)),
      norm(other.norm),
      order(other.order),
      childrenIDs(std::move(other.childrenIDs)),
      covered(std::move(other.covered)),
      halfspaces(std::move(other.halfspaces)) {}

    QNode& operator=(QNode&& other) noexcept {
        if (this != &other) {
            nodeID = other.nodeID;
            parentID = other.parentID;
            mbr = std::move(other.mbr);
            norm = other.norm;
            order = other.order;
            childrenIDs = std::move(other.childrenIDs);
            covered = std::move(other.covered);
            halfspaces = std::move(other.halfspaces);
        }
        return *this;
    }


    // Metodo per ottenere l'ID del nodo
    int getNodeID() const { return nodeID; }
    int getParentID() const { return parentID; }

    // Check if the node is the root
    bool isRoot() const;

    // Check if the node is a leaf
    bool isLeaf() const;

    // Get the order of the node
    void setOrder() const;

    // Get the covered halfspaces
    std::vector<HalfSpace> getCovered() const;

    // Insert halfspaces into the node
    void insertHalfspaces(const std::array<std::vector<std::vector<double>>, 2>& masks, const std::vector<HalfSpace>& halfspaces);

    // Getters and setters
    const std::vector<std::array<double, 2>>& getMBR() const { return mbr; }
    bool isNorm() const { return norm; }
    void setNorm(bool norm) { this->norm = norm; }
    size_t getOrder() const { return order; }

    // Get the children IDs instead of actual nodes
    const std::vector<int>& getChildrenIDs() const { return childrenIDs; }
    void addChildID(int childID) { childrenIDs.push_back(childID); }

    const std::vector<HalfSpace>& getHalfspaces() const { return halfspaces; }
    void setHalfspaces(const std::vector<HalfSpace>& halfspaces) { this->halfspaces = halfspaces; }
    void clearHalfspaces();

    // Serializes the QNode to disk
    void saveToDisk(const std::string& filePath);

    // Loads the QNode from disk
    void loadFromDisk(const std::string& filePath);

private:
    int nodeID{};  // ID univoco del nodo
    int parentID{};             // ID del nodo genitore

    std::vector<std::array<double, 2>> mbr;  // Minimum bounding region
    bool norm{};  // Normalization flag
    mutable size_t order{};  // Order of the node
    std::vector<int> childrenIDs;
    std::vector<HalfSpace> covered;  // Covered halfspaces
    std::vector<HalfSpace> halfspaces;  // Halfspaces in the node
};

#endif //QNODE_H
