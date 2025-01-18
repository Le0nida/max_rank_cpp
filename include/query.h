#ifndef QUERY_H
#define QUERY_H

#include <vector>
#include "geom.h"

/**
 * \brief Finds all points that strictly dominate point p in each dimension.
 * \param data A list of points.
 * \param p    The reference point.
 * \return A vector of dominating points.
 */
std::vector<Point> getdominators(const std::vector<Point>& data,
                                 const Point& p);

/**
 * \brief Finds all points that are dominated by p in each dimension.
 * \param data A list of points.
 * \param p    The reference point.
 * \return A vector of dominee points.
 */
std::vector<Point> getdominees(const std::vector<Point>& data,
                               const Point& p);

/**
 * \brief Finds points that are neither dominated by p nor dominating p.
 * \param data A list of points.
 * \param p    The reference point.
 * \return A vector of incomparable points.
 */
std::vector<Point> getincomparables(const std::vector<Point>& data,
                                    const Point& p);

/**
 * \brief Returns a set of points that form the skyline (not dominated by any other).
 * \param data The input set of points.
 * \return A vector of points in the skyline.
 */
std::vector<Point> getskyline(const std::vector<Point>& data);

#endif // QUERY_H
