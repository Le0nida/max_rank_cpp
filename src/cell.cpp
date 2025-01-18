#include "cell.h"
#include "halfspace.h"
#include "qtree.h"
#include <cmath>
#include <bitset>
#include <iostream>
#include <limits>
#include <iomanip>
#include <random>
#include <algorithm>
#include <utility>
#include <src/Highs.h>

/**
 * \class LinprogResult
 * \brief Internal structure to store results from HiGHS linear programming.
 */
struct LinprogResult {
    double* x;          ///< Pointer to the solution array
    double fun;         ///< Objective value
    int status;         ///< Solver status code
    char message[128];  ///< Solver status message
};

/**
 * \brief Frees a LinprogResult.
 */
static void free_linprog_result(LinprogResult* result) {
    delete[] result->x;
    delete result;
}

/// -------------------------------------------------
///                   Cell methods
/// -------------------------------------------------

Cell::Cell(const int order,
           std::string  mask,
           const std::vector<long>& covered,
           const std::vector<long>& halfspaces,
           const std::vector<std::array<double, 2>>& leaf_mbr,
           Point  feasible_pnt)
    : order(order),
      mask(std::move(mask)),
      covered(covered),
      halfspaces(halfspaces),
      leaf_mbr(leaf_mbr),
      feasible_pnt(std::move(feasible_pnt))
{
}

bool Cell::issingular() const {
    // Checks if all halfspaces in 'covered' are SINGULAR
    return std::all_of(covered.begin(), covered.end(),
        [](long id) {
            return halfspaceCache->get(id)->arr == Arrangement::SINGULAR;
        }
    );
}

/// -------------------------------------------------
///           Hamming Strings Generation
/// -------------------------------------------------

std::vector<std::string> genhammingstrings(int strlen, int weight) {
    // Limit bitstring length
    if (strlen > halfspacesLengthLimit) strlen = halfspacesLengthLimit;

    // Decide top-down or bottom-up
    bool botup;
    if (weight > std::ceil(strlen / 2.0)) {
        weight = strlen - weight;
        botup = false;
    } else {
        botup = true;
    }

    std::vector<int> decstr;
    decstr.reserve(1 << strlen);

    // Special cases for weight 0 or 1
    if (weight == 0) {
        decstr.push_back(0);
    }
    else if (weight == 1) {
        decstr.reserve(strlen);
        for (int b = 0; b < strlen; ++b) {
            decstr.push_back(1 << b);
        }
    }
    else {
        int halfmax = (1 << (strlen - 1)) - 1;
        int curr_weight = 2;

        // Start with bit patterns of weight=2 in range
        for (int b = 1; b < strlen; ++b) {
            decstr.push_back((1 << b) + 1);
        }
        std::vector<int> bases = decstr;
        bases.erase(std::remove_if(bases.begin(),
                                   bases.end(),
                                   [halfmax](int x) { return x > halfmax; }),
                    bases.end());

        // Expand until we reach the desired weight
        while (true) {
            // Expand horizontally
            while (!bases.empty()) {
                std::vector<int> shifts;
                for (int base : bases) {
                    int shifted = base << 1;
                    if (shifted <= halfmax) {
                        shifts.push_back(shifted);
                    }
                }
                decstr.insert(decstr.end(), shifts.begin(), shifts.end());
                bases = std::move(shifts);
            }

            // Increase weight if needed
            if (curr_weight < weight) {
                std::vector<int> new_bases;
                for (int dec : decstr) {
                    int new_dec = (dec << 1) + 1;
                    if (new_dec <= ((1 << strlen) - 1)) {
                        new_bases.push_back(new_dec);
                    }
                }
                decstr.insert(decstr.end(), new_bases.begin(), new_bases.end());
                bases = std::move(new_bases);
                curr_weight++;
            } else {
                break;
            }
        }
    }

    // Convert to actual bitstrings
    std::vector<std::string> hamming_strings;
    hamming_strings.reserve(decstr.size());
    for (int dec : decstr) {
        hamming_strings.emplace_back(
            std::bitset<64>(dec).to_string().substr(64 - strlen)
        );
    }

    // If top-down, invert bits
    if (!botup) {
        int decmax = (1 << strlen) - 1;
        for (auto& hs : hamming_strings) {
            int dec_val = static_cast<int>(
                std::bitset<64>(hs).to_ulong()
            );
            hs = std::bitset<64>(decmax - dec_val)
                    .to_string()
                    .substr(64 - strlen);
        }
    }

    return hamming_strings;
}

/// -------------------------------------------------
///          HiGHS Utility Functions
/// -------------------------------------------------

/**
 * \brief Calls the HiGHS solver with raw arrays in CSR format.
 */
static LinprogResult* linprog_highs(const double* c,
                                    const double* A_ub,
                                    const double* b_ub,
                                    const double* bounds,
                                    int num_vars,
                                    int num_constraints)
{
    auto* result = new LinprogResult;
    Highs highs;

    // Disable console output
    HighsOptions options;
    options.output_flag = false;
    highs.passOptions(options);

    // Problem dimensions
    const int num_col = num_vars;
    const int num_row = num_constraints;

    // Objective function
    std::vector<double> col_cost(c, c + num_col);

    // Build CSR data structures
    std::vector<int> A_start(num_row + 1);
    std::vector<int> A_index;
    A_index.reserve(num_row * num_col);
    std::vector<double> A_value;
    A_value.reserve(num_row * num_col);

    int count = 0;
    for (int i = 0; i < num_row; ++i) {
        A_start[i] = count;
        for (int j = 0; j < num_col; ++j) {
            double value = A_ub[i * num_col + j];
            if (std::fabs(value) > 1e-15) {
                A_index.push_back(j);
                A_value.push_back(value);
                count++;
            }
        }
    }
    A_start[num_row] = count;

    // Upper/lower constraints
    std::vector<double> row_upper(b_ub, b_ub + num_row);
    std::vector<double> row_lower(num_row, -kHighsInf);

    // Variable bounds
    std::vector<double> col_lower(num_col);
    std::vector<double> col_upper(num_col);

    for (int i = 0; i < num_col; ++i) {
        col_lower[i] = bounds[2 * i];
        col_upper[i] = bounds[2 * i + 1];
    }

    // Add columns and rows
    highs.addCols(num_col, col_cost.data(), col_lower.data(), col_upper.data(),
                  0, nullptr, nullptr, nullptr);
    highs.addRows(num_row, row_lower.data(), row_upper.data(),
                  (int)A_value.size(), A_start.data(), A_index.data(), A_value.data());

    // Run solver
    highs.run();

    // Store solution
    HighsSolution solution = highs.getSolution();
    HighsModelStatus model_status = highs.getModelStatus();

    result->x = new double[num_col];
    for (int i = 0; i < num_col; ++i) {
        result->x[i] = solution.col_value[i];
    }
    result->fun = highs.getObjectiveValue();
    result->status = static_cast<int>(model_status);
    strncpy(result->message, highs.modelStatusToString(model_status).c_str(), 128);

    return result;
}

/**
 * \brief Wraps linprog_highs to handle std::vector-based inputs, returning a tuple.
 */
std::tuple<std::vector<double>, double, int, std::string>
linprog_highs(const std::vector<double>& c,
              const std::vector<std::vector<double>>& A_ub,
              const std::vector<double>& b_ub,
              const std::vector<std::pair<double, double>>& bounds)
{
    int num_vars = (int)c.size();
    int num_constraints = (int)b_ub.size();

    // Flatten inputs
    std::vector<double> c_array(c.begin(), c.end());

    std::vector<double> A_ub_array;
    A_ub_array.reserve((size_t)num_constraints * num_vars);
    for (const auto& row : A_ub) {
        A_ub_array.insert(A_ub_array.end(), row.begin(), row.end());
    }

    std::vector<double> b_ub_array(b_ub.begin(), b_ub.end());

    std::vector<double> bounds_array;
    bounds_array.reserve((size_t)2 * num_vars);
    for (const auto& bound : bounds) {
        bounds_array.push_back(bound.first);
        bounds_array.push_back(bound.second);
    }

    // Solve
    LinprogResult* raw_result = linprog_highs(c_array.data(),
                                              A_ub_array.data(),
                                              b_ub_array.data(),
                                              bounds_array.data(),
                                              num_vars,
                                              num_constraints);

    // Extract result
    std::vector<double> solution(raw_result->x, raw_result->x + num_vars);
    double fun = raw_result->fun;
    int status = raw_result->status;
    std::string message(raw_result->message);

    free_linprog_result(raw_result);
    return {solution, fun, status, message};
}

/// -------------------------------------------------
///       Search minimal cells using LP
/// -------------------------------------------------

std::vector<Cell> searchmincells_lp(const QNode& leaf,
                                    const std::vector<std::string>& hamstrings)
{
    std::vector<Cell> cells;
    cells.reserve(hamstrings.size());

    int dims = (int)leaf.mbr.size();
    std::vector<long> leaf_covered = leaf.getCovered();
    auto halfspaces = leaf.halfspaces;

    // If no halfspaces, build a trivial cell from MBR center
    if (halfspaces.empty()) {
        std::vector<std::array<double, 2>> mbr = leaf.mbr;
        std::vector<double> center(dims);
        for (int i = 0; i < dims; ++i) {
            center[i] = 0.5 * (mbr[i][0] + mbr[i][1]);
        }
        Point feasible_pnt(center);
        return { Cell(0, "", leaf_covered, {}, mbr, feasible_pnt) };
    }

    // Build objective with slack variable
    std::vector<double> c(dims + 1, 0.0);
    c[dims] = -1.0; // Minimizing negative slack => maximizing slack

    // Prepare A_ub
    std::vector<std::vector<double>> A_ub(halfspaces.size() + 1,
                                          std::vector<double>(dims + 1, 1.0));
    A_ub[halfspaces.size()][dims] = 0.0;

    std::vector<double> b_ub(halfspaces.size() + 1, 1.0);

    // Bounds for each variable + slack
    std::vector<std::pair<double,double>> bounds(dims + 1);
    for (int d = 0; d < dims; ++d) {
        bounds[d] = {leaf.mbr[d][0], leaf.mbr[d][1]};
    }
    bounds[dims] = {0.0, std::numeric_limits<double>::infinity()};

    // Try each Hamming string
    for (const auto& hamstr : hamstrings) {
        // Set constraints according to the bitstring
        for (int b = 0; b < (int)hamstr.size(); ++b) {
            auto hs = halfspaceCache->get(halfspaces[b]);
            if (!hs) continue; // Just in case

            if (hamstr[b] == '0') {
                for (int i = 0; i < dims; ++i) {
                    A_ub[b][i] = -hs->coeff[i];
                }
                b_ub[b] = -hs->known;
            } else {
                for (int i = 0; i < dims; ++i) {
                    A_ub[b][i] = hs->coeff[i];
                }
                b_ub[b] = hs->known;
            }
        }

        // Solve the LP
        auto [solution, fun, status, message] = linprog_highs(c, A_ub, b_ub, bounds);

        // If feasible, build a Cell
        if (status == (int)HighsModelStatus::kOptimal) {
            Point feasible_pnt(std::vector<double>(solution.begin(),
                                                   solution.end() - 1));
            cells.emplace_back(0, hamstr, leaf_covered,
                               halfspaces, leaf.mbr, feasible_pnt);
            // Break after first feasible
            break;
        }
    }

    return cells;
}
