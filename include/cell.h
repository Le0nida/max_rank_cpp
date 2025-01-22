#ifndef CELL_H
#define CELL_H

#include <algorithm>
#include <utility>
#include <vector>
#include <string>
#include <array>
#include <memory>
#include "geom.h"
#include "halfspace.h"
#include "qnode.h"

extern int halfspacesLengthLimit;
extern int maxNoBinStringToCheck;

/**
 * \class Cell
 * \brief Represents a minimal region in the halfspace expansion process.
 */
class Cell {
public:
    /**
     * \brief Constructor
     * \param order        The "order" of this cell (sum of halfspaces).
     * \param mask         A binary mask representing halfspace configuration.
     * \param covered      List of halfspace IDs fully covering this cell.
     * \param halfspaces   List of overlapping halfspaces.
     * \param leaf_mbr     The [min,max] bounds in each dimension for this cell.
     * \param feasible_pnt A point guaranteed to lie in this cell.
     */
    Cell(int order,
         std::string  mask,
         const std::vector<long>& covered,
         const std::vector<long>& halfspaces,
         const std::vector<std::array<float, 2>>& leaf_mbr,
         Point  feasible_pnt);

    /**
     * \brief Checks if all halfspaces covering this cell are marked as SINGULAR.
     * \return True if every halfspace in 'covered' has Arrangement::SINGULAR.
     */
    [[nodiscard]] bool issingular() const;

    int order;                                   ///< Summed order for this cell
    std::string mask;                            ///< Binary string representation
    std::vector<long> covered;                   ///< Fully covered halfspaces
    std::vector<long> halfspaces;                ///< Overlapping halfspaces
    std::vector<std::array<float, 2>> leaf_mbr; ///< Local bounding region
    Point feasible_pnt;                          ///< A feasible point in this region
};

/**
 * \struct Interval
 * \brief Used in 2D expansions, representing a range on the x-axis with coverage info.
 */
struct Interval {

    std::shared_ptr<HalfLine> halfline;             ///< A reference to the halfline generating this interval
    std::pair<double, double> range;                ///< Range of x-coordinates [range.first, range.second]
    bool coversleft;                                ///< If true, this halfline "exits" at range.second; otherwise it "enters"
    int order;                                      ///< Number of halflines covering this interval
    std::vector<std::shared_ptr<HalfLine>> covered; ///< Collection of halflines covering this interval

    /**
     * \brief Constructor for Interval
     */
    Interval(std::shared_ptr<HalfLine> hl,
             std::pair<double, double> r,
             const bool c)
        : halfline(std::move(hl)), range(std::move(r)), coversleft(c), order(0) {}

    /**
     * \brief Checks if all halflines covering this interval are SINGULAR.
     * \return True if all are Arrangement::SINGULAR.
     */
    [[nodiscard]] bool issingular() const {
        return std::all_of(covered.begin(), covered.end(),
            [](const std::shared_ptr<HalfLine>& hl) {
                return hl->arr == Arrangement::SINGULAR;
            }
        );
    }
};

/**
 * \brief Generates Hamming strings of a specific weight.
 * \param strlen Length of the bitstring.
 * \param weight Desired Hamming weight.
 * \return A list of bitstrings.
 */
std::vector<std::string> genhammingstrings(int strlen, int weight);

/**
 * \brief Searches for minimal cells using linear programming.
 * \param leaf        A reference to a QNode (leaf) with bounding MBR and halfspaces.
 * \param hamstrings  A list of precomputed Hamming strings.
 * \return A list of Cell objects that pass the feasibility check.
 */
std::vector<Cell> searchmincells_lp(const QNode& leaf,
                                    const std::vector<std::string>& hamstrings);

#endif // CELL_H
