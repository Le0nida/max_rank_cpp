//
// Created by leona on 01/08/2024.
//

#ifndef QNODE_H
#define QNODE_H

#include <utility>
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


    QNode(QNode&& other) noexcept {}

    QNode& operator=(QNode&& other) noexcept {return *this;}


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
    std::vector<long int> getCovered() const;

    // Insert halfspaces into the node
    void insertHalfspaces(const std::array<std::vector<std::vector<double>>, 2>& masks, const std::vector<long int>& halfspaces);

    // Getters and setters
    const std::vector<std::array<double, 2>>& getMBR() const { return mbr; }
    bool isNorm() const { return norm; }
    void setNorm(bool norm) { this->norm = norm; }
    size_t getOrder() const { return order; }

    // Get the children IDs instead of actual nodes
    const std::vector<long int>& getChildrenIDs() const { return childrenIDs; }
    void addChildID(int childID) { childrenIDs.push_back(childID); }

    const std::vector<long int>& getHalfspaces() const { return halfspaces; }
    void setHalfspaces(std::vector<long int> halfspaces) { this->halfspaces = std::move(halfspaces); };
    void clearHalfspaces();


    void saveToDisk(const std::string& filePath);           // Serializes the QNode to disk
    void loadFromDisk(const std::string& filePath);         // Loads the QNode from disk

private:
    long int nodeID;                        // ID univoco del nodo
    long int parentID;                      // ID del nodo genitore

    std::vector<std::array<double, 2>> mbr; // Minimum bounding region
    bool norm;                              // Normalization flag
    mutable size_t order{};                 // Order of the node

    std::vector<long int> childrenIDs;
    std::vector<long int> covered;          // Covered halfspaces
    std::vector<long int> halfspaces;       // Halfspaces in the node
};

#endif //QNODE_H
