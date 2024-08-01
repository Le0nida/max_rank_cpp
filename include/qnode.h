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

    bool isRoot() const;
    bool isLeaf() const;
    int getOrder();
    std::vector<Halfspace> getCovered();
    void insertHalfspaces(const std::vector<std::array<std::vector<double>, 2>>& masks, const std::vector<Halfspace>& halfspaces);

    const std::vector<std::array<double, 2>>& getMBR() const { return mbr; }
    bool isNorm() const { return norm; }
    void setNorm(bool norm) { this->norm = norm; }
    const std::vector<std::unique_ptr<QNode>>& getChildren() const { return children; }
    void addChild(std::unique_ptr<QNode> child) { children.push_back(std::move(child)); }
    const std::vector<Halfspace>& getHalfspaces() const { return halfspaces; }
    void setHalfspaces(const std::vector<Halfspace>& halfspaces) { this->halfspaces = halfspaces; }
    void clearHalfspaces() { halfspaces.clear(); }

private:
    std::vector<std::array<double, 2>> mbr;
    bool norm;
    int order;
    QNode* parent;
    std::vector<std::unique_ptr<QNode>> children;
    std::vector<Halfspace> covered;
    std::vector<Halfspace> halfspaces;
};

#endif //QNODE_H
