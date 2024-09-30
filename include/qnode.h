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

extern int globalNodeID;

enum class PositionHS { BELOW, ABOVE, OVERLAPPED };

class QNode {
public:
    QNode(int parentID, const std::vector<std::array<double, 2>>& mbr);
    QNode();
    ~QNode() = default;
    explicit QNode(int id) : nodeID(id) {}

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
    //void insertHalfspaces(const std::array<std::vector<std::vector<double>>, 2>& masks, const std::vector<long int>& halfspaces);
    void insertHalfspaces(const std::vector<long int>& halfspaces);

    PositionHS MbrVersusHalfSpace(const std::vector<double>& hs_coeff, double hs_known, const std::vector<std::array<double, 2>>& mbr);

    // Getters and setters
    const std::vector<std::array<double, 2>>& getMBR() const { return mbr; }
    bool isNorm() const { return norm; }
    void setNorm(bool norm) { this->norm = norm; }
    void setLeaf(bool leaf) { this->leaf = leaf; }
    size_t getOrder() const { return order; }

    // Get the children IDs instead of actual nodes
    const std::vector<long int>& getChildrenIDs() const { return childrenIDs; }
    void addChildID(int childID) { childrenIDs.push_back(childID); }

    const std::vector<long int>& getHalfspaces() const { return halfspaces; }
    void setHalfspaces(std::vector<long int> halfspaces) { this->halfspaces = std::move(halfspaces); };
    void clearHalfspaces();


    size_t saveToDisk(std::fstream& outStream);           // Serializes the QNode to disk
    void loadFromDisk(std::fstream& inStream);         // Loads the QNode from disk

    void splitNode();
    std::vector<std::vector<std::array<double, 2>>> genSubdivisions();
    bool checkNodeValidity();

private:
    long int nodeID;                        // ID univoco del nodo
    long int parentID;                      // ID del nodo genitore

    std::vector<std::array<double, 2>> mbr; // Minimum bounding region
    bool norm;                              // Normalization flag
    bool leaf{};                              // Leaf flag
    mutable size_t order{};                 // Order of the node

    std::vector<long int> childrenIDs;
    std::vector<long int> covered;          // Covered halfspaces
    std::vector<long int> halfspaces;       // Halfspaces in the node

    size_t estimatedSize() const;
};

#endif //QNODE_H
