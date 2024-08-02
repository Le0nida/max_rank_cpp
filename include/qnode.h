//
// Created by leona on 01/08/2024.
//

#ifndef QNODE_H
#define QNODE_H

#include <vector>
#include <array>
#include <memory>

class Halfspace;

class QNode {
public:
    QNode(QNode* parent, const std::vector<std::array<double, 2>>& mbr);
    ~QNode() = default;

    // Check if the node is the root
    bool isRoot() const;

    // Check if the node is a leaf
    bool isLeaf() const;

    // Get the order of the node
    size_t getOrder() const;

    // Get the covered halfspaces
    std::vector<Halfspace> getCovered() const;

    // Insert halfspaces into the node
    void insertHalfspaces(const std::vector<std::array<std::vector<double>, 2>>& masks, const std::vector<Halfspace>& halfspaces);

    // Getters and setters
    const std::vector<std::array<double, 2>>& getMBR() const { return mbr; }
    bool isNorm() const { return norm; }
    void setNorm(bool norm) { this->norm = norm; }
    const std::vector<std::unique_ptr<QNode>>& getChildren() const { return children; }
    void addChild(std::unique_ptr<QNode> child) { children.push_back(std::move(child)); }
    const std::vector<Halfspace>& getHalfspaces() const { return halfspaces; }
    void setHalfspaces(const std::vector<Halfspace>& halfspaces) { this->halfspaces = halfspaces; }
    void clearHalfspaces() { halfspaces.clear(); }

private:
    std::vector<std::array<double, 2>> mbr;  // Minimum bounding region
    bool norm;  // Normalization flag
    mutable size_t order;  // Order of the node
    QNode* parent;  // Parent node
    std::vector<std::unique_ptr<QNode>> children;  // Children nodes
    std::vector<Halfspace> covered;  // Covered halfspaces
    std::vector<Halfspace> halfspaces;  // Halfspaces in the node
};

#endif //QNODE_H
