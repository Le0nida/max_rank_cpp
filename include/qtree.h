#ifndef QTREE_H
#define QTREE_H

#include <vector>
#include <array>
#include "qnode.h"
#include <algorithm>
#include <future>
#include <iostream>
#include <queue>
#include <thread>

/**
 * \class QTree
 * \brief Manages an N-dimensional tree structure (like a quadtree or octree)
 *        for storing and subdividing halfspaces.
 */
class QTree {
public:

    int dims;       ///< Number of dimensions in the reduced query space
    int maxhsnode;  ///< Max halfspaces per node before triggering a split
    int maxLevel;   ///< Maximum depth allowed in this tree

    QNode* root;                     ///< Root node

    std::vector<QNode*> macroRoots;  ///< Collection of macro-root nodes (one per sub-MBR)
    /**
     * \brief Precomputed subdivisions of the unit hypercube:
     *        each sub-MBR is a vector of [min,max] pairs in float.
     */
    std::vector< std::vector<std::array<float,2>> > precomputedSubMBRs;

    /**
     * \brief Constructor
     * \param dims       Number of dimensions.
     * \param maxhsnode  Max halfspaces per node (splitting threshold).
     * \param maxLevel   Maximum allowed tree depth.
     */
    QTree(int dims, int maxhsnode, int maxLevel);

    /**
     * \brief Destructor. Destroys the root and all macro-roots.
     */
    ~QTree();

    /**
     * \brief Builds an empty root node covering [0,1]^dims.
     * \return A pointer to the newly created QNode.
     */
    QNode* createroot();

    /**
     * \brief Destroys both the classical root and all macro-root subtrees.
     */
    void destroyAllNodes();

    /**
     * \brief Inserts halfspaces in parallel, splitting them among the precomputed sub-MBRs.
     *        If a sub-MBR is fully covered (no partial overlap), no macro-root is created.
     * \param halfspaces Vector of halfspace IDs to distribute and insert.
     */
    void inserthalfspacesMacroSplit(const std::vector<long int>& halfspaces);

    /**
     * \brief Returns all leaves from all macro-root subtrees.
     * \return Vector of leaf pointers.
     */
    [[nodiscard]] std::vector<QNode*> getAllLeaves() const;

    /**
     * \brief Recomputes the order (rank) of each node in the root and all macro-roots.
     */
    void updateAllOrders();

    /**
     * \brief Builds a new subtree for the given sub-MBR and halfspace set.
     * \param subMBR  The bounding region for this subtree (float).
     * \param subHS   Vector of halfspace IDs to insert incrementally.
     * \return The newly created QNode pointer as a subtree root.
     */
    [[nodiscard]] QNode* buildSubtree(const std::vector<std::array<float,2>>& subMBR,
                                      const std::vector<long>& subHS) const;

    /**
     * \brief Precomputes 2^dims sub-MBRs covering [0,1]^dims, stored as float.
     * \param globalMBR The bounding region, typically the entire unit hypercube (float).
     * \return A list of subdivided MBRs in float.
     */
    [[nodiscard]] std::vector< std::vector<std::array<float,2>> >
    macroSplitMBR(const std::vector<std::array<float,2>>& globalMBR) const;
};

#endif // QTREE_H
