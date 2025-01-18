#ifndef HALFSPACE_H
#define HALFSPACE_H

#include <memory>
#include <unordered_map>
#include <vector>
#include <fstream>
#include "geom.h"

/**
 * \enum Position
 * \brief Spatial position of a point relative to a halfspace.
 */
enum Position {
    POS_IN  = 1,  ///< Point is inside or below the halfspace
    POS_OUT = -1, ///< Point is outside or above
    POS_ON  = 0   ///< Point is exactly on the boundary
};

/**
 * \enum Arrangement
 * \brief Describes how a halfspace or halfline is marked during expansions.
 */
enum class Arrangement {
    SINGULAR   = 0, ///< Marked as finalized/expanded
    AUGMENTED  = 1  ///< Still subject to further expansion
};

/**
 * \class HalfLine
 * \brief In 2D, represents a line y = m*x + q from a pivot point.
 */
class HalfLine {
public:
    /**
     * \brief Constructs a halfline from a given point (2D).
     * \param pnt The reference Point with x,y coordinates.
     */
    explicit HalfLine(const Point& pnt);

    /**
     * \brief Computes y = m*x + q
     * \param x The x-coordinate
     * \return The corresponding y-value on this halfline.
     */
    [[nodiscard]] double get_y(double x) const;

    Point pnt;           ///< Reference point used to build the line
    double m;            ///< Slope
    double q;            ///< Intercept
    Arrangement arr;     ///< Indicates if it's SINGULAR or AUGMENTED
    int dims;            ///< Typically 2 for 2D lines
};

/**
 * \class HalfSpace
 * \brief Represents a linear halfspace in N dimensions.
 */
class HalfSpace {
public:
    /**
     * \brief Constructs a halfspace from an ID, coefficients, and known term.
     * \param pntID Identifier for the halfspace (often the point ID).
     * \param coeff The vector of coefficients in the linear inequality.
     * \param known The RHS constant.
     */
    HalfSpace(long int pntID,
              const std::vector<double>& coeff,
              double known);

    /**
     * \brief Default constructor (invalid halfspace).
     */
    HalfSpace();

    long int pntID;            ///< Associated point ID
    std::vector<double> coeff; ///< Coefficients in each dimension
    double known;              ///< RHS constant
    Arrangement arr;           ///< Mark as AUGMENTED or SINGULAR
    int dims;                  ///< Number of dimensions

    /**
     * \brief Equality operator
     */
    bool operator==(const HalfSpace& other) const;
};

/**
 * \brief Finds the intersection between two 2D halflines r and s.
 * \return Point with intersection coords, or (inf, inf) if parallel.
 */
Point find_halflines_intersection(const HalfLine& r, const HalfLine& s);

/**
 * \brief Checks the position of a Point relative to a HalfSpace.
 * \param point The point to test.
 * \param halfspace The halfspace.
 * \return POS_IN, POS_OUT, or POS_ON.
 */
Position find_pointhalfspace_position(const Point& point, const HalfSpace& halfspace);

/**
 * \brief Generates halfspaces associated with point p for each record in 'records'.
 * \param p       Reference point.
 * \param records A list of points to convert.
 * \return A vector of halfspace IDs.
 */
std::vector<long> genhalfspaces(const Point& p,
                                const std::vector<Point>& records);

/**
 * \class HalfSpaceCache
 * \brief A global cache storing HalfSpace objects by their IDs.
 */
class HalfSpaceCache {
public:
    /**
     * \brief Constructor preallocating a given number of entries.
     * \param cacheSize Expected maximum number of halfspaces.
     */
    explicit HalfSpaceCache(size_t cacheSize) {
        cache.reserve(cacheSize);
    }

    /**
     * \brief Inserts or updates a halfspace with the given ID.
     * \param id         Halfspace ID
     * \param halfspace  Shared pointer to the HalfSpace.
     */
    void insert(long id, std::shared_ptr<HalfSpace> halfspace) {
        cache[id] = std::move(halfspace);
    }

    /**
     * \brief Retrieves a shared pointer to the HalfSpace with the given ID.
     * \param id The halfspace ID.
     * \return A valid shared_ptr if found, or null if not found.
     */
    std::shared_ptr<HalfSpace> get(long id) const {
        auto it = cache.find(id);
        if (it != cache.end()) {
            return it->second;
        }
        return nullptr;
    }

    /**
     * \brief Checks if a halfspace with the given ID is present in the cache.
     * \param id The halfspace ID.
     * \return True if found, otherwise false.
     */
    bool contains(long id) const {
        return (cache.find(id) != cache.end());
    }

private:
    /**
     * \brief Internal storage of <halfspaceID, pointer> pairs.
     */
    std::unordered_map<long, std::shared_ptr<HalfSpace>> cache;
};

/**
 * \struct PointHash
 * \brief Hash functor for using Point as a key in std::unordered_map.
 */
struct PointHash {
    std::size_t operator()(const Point& p) const {
        std::size_t seed = 0;
        for (double coord : p.coord) {
            // Combine using a standard "hash mix" formula
            seed ^= std::hash<double>{}(coord) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

/**
 * \brief Global pointer to the halfspace cache. Initialize via initializeCache().
 */
extern HalfSpaceCache* halfspaceCache;

/**
 * \brief Initializes halfspaceCache with a fixed size.
 * \param cacheSize The maximum number of halfspaces to store.
 */
static void initializeCache(const size_t cacheSize) {
    if (!halfspaceCache) {
        halfspaceCache = new HalfSpaceCache(cacheSize);
    }
}

/**
 * \brief Maps a Point to its halfspace ID if already converted.
 */
extern std::unordered_map<Point, long, PointHash> pointToHalfSpaceCache;

#endif // HALFSPACE_H
