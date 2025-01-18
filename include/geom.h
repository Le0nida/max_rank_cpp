#ifndef GEOM_H
#define GEOM_H

#include <vector>
#include <fstream>
#include <string>

/**
 * \class Point
 * \brief A generic point in N-dimensional space.
 */
class Point {
public:
    /**
     * \brief Constructs a point from coordinates and an optional ID.
     * \param coord Vector of coordinates.
     * \param id    Optional integer identifier.
     */
    explicit Point(const std::vector<double>& coord, int id = -1);

    int id;                      ///< Unique ID (if any)
    std::vector<double> coord;   ///< The point coordinates
    int dims;                    ///< Number of dimensions in coord

    /**
     * \brief Equality operator checks if coordinates match exactly.
     */
    bool operator==(const Point& other) const {
        return (coord == other.coord);
    }
};

/**
 * \brief Generates masks for geometric operations (demo usage).
 * \param dims Number of dimensions.
 * \return Pair of mask vectors.
 */
std::pair<std::vector<std::vector<double>>,
          std::vector<std::vector<double>>> genmasks(int dims);

#endif // GEOM_H
