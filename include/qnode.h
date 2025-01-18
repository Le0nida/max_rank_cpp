#ifndef QNODE_H
#define QNODE_H

#include <vector>
#include <array>
#include "halfspace.h"

extern int numOfSubdivisions;    ///< Global number of partitions used in node splitting (2^dims)

/**
 * \enum PositionHS
 * \brief Relative position of a node's MBR with respect to a halfspace.
 */
enum class PositionHS { BELOW, ABOVE, OVERLAPPED };

class QTree;

/**
 * \class QNode
 * \brief Represents a node in the QTree. Each node can store halfspaces
 *        and subdivide into children if necessary.
 *
 * All members are public for simplicity. MBR covers [min,max] in each dimension.
 */
class QNode {
public:

    QTree* owner;                 ///< Pointer to the owner QTree
    QNode* parent;                ///< Pointer to the parent QNode (null if root)
    std::vector<QNode*> children; ///< Child nodes (size = 2^dims). Some may be null if not valid.
    std::vector<long> covered;    ///< Halfspaces that fully cover this node (no further checks required)
    std::vector<long> halfspaces; ///< Halfspaces that overlap this node and need potential splitting


    int leafIndex;                          ///< Index used if needed (e.g., in a list of leaves)
    std::vector<std::array<double,2>> mbr;  ///< [min,max] bounding region for each dimension

    bool norm;     ///< True if this node is valid
    bool leaf;     ///< True if this node is a leaf (no children)
    size_t order;  ///< Accumulated order (sum of covered halfspaces of parents)
    int level;     ///< Depth level in the tree (0 = root)

    /**
     * \brief Constructor
     * \param owner  Owning QTree.
     * \param parent Parent QNode (or null if root).
     * \param mbr    Bounding rectangle for each dimension.
     * \param level  Depth level in the QTree.
     */
    QNode(QTree* owner,
          QNode* parent,
          const std::vector<std::array<double,2>>& mbr,
          int level);

    ///< Delete copy semantics to avoid accidental duplication
    QNode(const QNode&) = delete;
    QNode& operator=(const QNode&) = delete;

    /**
     * \brief Sets the leaf status of this node.
     */
    void setLeaf(bool lf);

    /**
     * \brief Inserts multiple halfspaces into this node, one by one.
     * \param new_halfspaces Vector of halfspace IDs.
     */
    void insertHalfspaces(const std::vector<long>& new_halfspaces);

    /**
     * \brief Inserts a single halfspace into this node, possibly triggering splitting.
     * \param hsID Identifier of the halfspace.
     */
    void insertHalfspace(long hsID);

    /**
     * \brief Splits this node into children if capacity is exceeded.
     */
    void splitNode();

    /**
     * \brief Checks if the node is valid in [0,1]^dims (sum of corners <= 1).
     * \return True if valid, false otherwise.
     */
    [[nodiscard]] bool checkNodeValidity() const;

    /**
     * \brief Clears the halfspaces vector in this node.
     */
    void clearHalfspaces();

    /**
     * \brief Determines how the node's MBR relates to the specified halfspace.
     * \param coeff Coefficients of the halfspace.
     * \param known Right-hand side (RHS) of the halfspace.
     * \return PositionHS::BELOW, ABOVE, or OVERLAPPED.
     */
    [[nodiscard]] PositionHS MbrVersusHalfSpace(const std::vector<double>& coeff, double known) const;

    /**
     * \brief Gathers 'covered' halfspaces from this node and all ancestor nodes.
     * \return A vector of halfspace IDs.
     */
    [[nodiscard]] std::vector<long> getCovered() const;

    /**
     * \brief Checks if this node is the root.
     * \return True if it has no parent.
     */
    [[nodiscard]] bool isRoot() const { return (parent == nullptr); }
};

#endif // QNODE_H
