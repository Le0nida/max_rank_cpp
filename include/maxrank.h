#ifndef MAXRANK_H
#define MAXRANK_H

#include <vector>
#include "geom.h"
#include "cell.h"
#include "query.h"
#include "halfspace.h"
#include "qtree.h"
#include "utils.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

/**
 * \brief Global configurable parameters for the MaxRank approach.
 */
extern int limitHamWeight;         ///< Max Hamming weight to consider
extern int maxLevelQTree;          ///< Maximum allowed QTree depth
extern int maxCapacityQNode;       ///< Maximum capacity of halfspaces in a QNode
extern int maxNoBinStringToCheck;  ///< Maximum number of binary strings to check

/**
 * \brief Main function for multi-dimensional MaxRank (d > 2).
 * \param data A set of points in the dataset.
 * \param p    The reference point for the MaxRank query.
 * \return Pair containing (MaxRank value, list of minimal Cells).
 *
 * \note This algorithm expands halfspaces around \p p, subdividing the space
 *       to find minimal cells that satisfy the ordering constraints.
 */
std::pair<int, std::vector<Cell>> aa_hd(const std::vector<Point>& data,
                                        const Point& p);

/**
 * \brief Special case for 2D MaxRank approach.
 * \param data A set of points in the dataset (2D).
 * \param p    The reference point in 2D.
 * \return Pair containing (MaxRank value, list of minimal Intervals).
 *
 * \note In 2D, the halfspace expansion is replaced by halfline expansions,
 *       producing intervals on the x-axis.
 */
std::pair<int, std::vector<Interval>> aa_2d(const std::vector<Point>& data,
                                            const Point& p);

#endif // MAXRANK_H
