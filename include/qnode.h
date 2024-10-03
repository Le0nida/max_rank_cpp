//
// Created by leona on 01/08/2024.
//

#ifndef QNODE_H
#define QNODE_H

#include <utility>
#include <vector>
#include <memory>
#include "geom.h"
#include "context.h"  // Include Context definition

class Context; // Forward declaration

enum class PositionHS { BELOW, ABOVE, OVERLAPPED };

class QNode {
public:
    QNode(Context& ctx, int parentID, const std::vector<std::array<double, 2>>& mbr);
    explicit QNode(Context& ctx);
    ~QNode() = default;

    QNode(const QNode&) = delete;
    QNode& operator=(const QNode&) = delete;

    QNode(QNode&& other) noexcept;
    QNode& operator=(QNode&& other) noexcept;

    int getNodeID() const;
    int getParentID() const;

    bool isRoot() const;
    bool isLeaf() const;

    void setOrder() const;

    std::vector<long int> getCovered() const;

    void insertHalfspaces(Context& ctx, const std::vector<long int>& halfspaces);

    PositionHS MbrVersusHalfSpace(const std::vector<double>& hs_coeff, double hs_known, const std::vector<std::array<double, 2>>& mbr);

    const std::vector<std::array<double, 2>>& getMBR() const;
    bool isNorm() const;
    void setNorm(bool norm);
    void setLeaf(bool leaf);
    size_t getOrder() const;

    const std::vector<long int>& getChildrenIDs() const;
    void addChildID(int childID);

    const std::vector<long int>& getHalfspaces() const;
    void setHalfspaces(std::vector<long int> halfspaces);
    void clearHalfspaces();

    void splitNode(Context& ctx);
    std::vector<std::vector<std::array<double, 2>>> genSubdivisions();
    bool checkNodeValidity();

    size_t saveToDisk(char* pFileData, size_t offset);
    void loadFromDisk(char* pFileData, size_t offset);
    size_t estimatedSize() const;

private:
    Context& ctx;
    long int nodeID;
    long int parentID;

    std::vector<std::array<double, 2>> mbr;
    bool norm;
    bool leaf;
    mutable size_t order;

    std::vector<long int> childrenIDs;
    std::vector<long int> covered;
    std::vector<long int> halfspaces;
};

#endif // QNODE_H
