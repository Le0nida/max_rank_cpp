//
// Created by leona on 06/08/2024.
//

#include "cell.h"
#include "halfspace.h"
#include "qtree.h"
#include <cmath>
#include <algorithm>
#include <random>
#include <Eigen/Dense>
#include <iostream>
#include <limits>
#include <iomanip>
#include <stdexcept>
#include <bitset>
#include <src/Highs.h>

Cell::Cell(int order, const std::string& mask, const std::vector<HalfSpace>& covered,
           const std::vector<HalfSpace>& halfspaces, const std::vector<std::array<double, 2>>& leaf_mbr,
           const Point& feasible_pnt)
    : order(order), mask(mask), covered(covered), halfspaces(halfspaces), leaf_mbr(leaf_mbr), feasible_pnt(feasible_pnt)
{
}

bool Cell::issingular() const
{
    return std::all_of(covered.begin(), covered.end(),
                       [](const HalfSpace& hs) { return hs.arr == Arrangement::SINGULAR; });
}

Interval::Interval(const HalfLine& halfline, const std::pair<double, double>& range, int coversleft)
    : halfline(halfline), range(range), coversleft(coversleft)
{
}

bool Interval::issingular() const
{
    return std::all_of(covered.begin(), covered.end(),
                       [](const HalfSpace& hl) { return hl.arr == Arrangement::SINGULAR; });
}

std::vector<std::string> genhammingstrings(int strlen, int weight)
{
    bool botup;
    if (weight > std::ceil(strlen / 2.0))
    {
        weight = strlen - weight;
        botup = false;
    }
    else
    {
        botup = true;
    }

    std::vector<int> decstr;
    if (weight == 0)
    {
        decstr = {0};
    }
    else if (weight == 1)
    {
        for (int b = 0; b < strlen; ++b)
        {
            decstr.push_back(1 << b);
        }
    }
    else
    {
        int halfmax = (1 << (strlen - 1)) - 1;
        int curr_weight = 2;
        for (int b = 1; b < strlen; ++b)
        {
            decstr.push_back((1 << b) + 1);
        }
        std::vector<int> bases = decstr;
        bases.erase(std::remove_if(bases.begin(), bases.end(), [halfmax](int x) { return x > halfmax; }), bases.end());

        while (true)
        {
            while (!bases.empty())
            {
                std::vector<int> shifts;
                for (int base : bases)
                {
                    int shifted = base << 1;
                    if (shifted <= halfmax)
                    {
                        shifts.push_back(shifted);
                    }
                }
                decstr.insert(decstr.end(), shifts.begin(), shifts.end());
                bases = shifts;
            }

            if (curr_weight < weight)
            {
                std::vector<int> new_bases;
                for (int dec : decstr)
                {
                    int new_dec = (dec << 1) + 1;
                    if (new_dec <= (1 << strlen) - 1)
                    {
                        new_bases.push_back(new_dec);
                    }
                }
                decstr.insert(decstr.end(), new_bases.begin(), new_bases.end());
                bases = new_bases;
                curr_weight++;
            }
            else
            {
                break;
            }
        }
    }

    std::vector<std::string> hamming_strings;
    for (int dec : decstr)
    {
        hamming_strings.push_back(std::bitset<32>(dec).to_string().substr(32 - strlen));
    }

    if (!botup)
    {
        int decmax = (1 << strlen) - 1;
        for (auto& hs : hamming_strings)
        {
            int dec_val = std::bitset<32>(hs).to_ulong();
            hs = std::bitset<32>(decmax - dec_val).to_string().substr(32 - strlen);
        }
    }

    return hamming_strings;
}

struct LinprogResult
{
    double* x;
    double fun;
    int status;
    char message[128];
};

void free_linprog_result(LinprogResult* result) {
    delete[] result->x;
    delete result;
}


LinprogResult* linprog_highs(const double* c, const double* A_ub, const double* b_ub,
                             const double* bounds, int num_vars, int num_constraints)
{
    auto* result = new LinprogResult;

    Highs highs;

    // Disabilita l'output a console di HiGHS
    HighsOptions options;
    options.output_flag = false;
    highs.passOptions(options);

    // Define the problem dimensions
    const int num_col = num_vars;
    const int num_row = num_constraints;

    // Objective function coefficients
    std::vector<double> col_cost(c, c + num_col);

    // Prepare A_ub in CSR format
    std::vector<int> A_start(num_row + 1);
    std::vector<int> A_index;
    std::vector<double> A_value;

    int count = 0;
    for (int i = 0; i < num_row; ++i)
    {
        A_start[i] = count;
        for (int j = 0; j < num_col; ++j)
        {
            double value = A_ub[i * num_col + j];
            if (value != 0.0)
            {
                A_index.push_back(j);
                A_value.push_back(value);
                count++;
            }
        }
    }
    A_start[num_row] = count;

    // Right-hand side values (upper bounds)
    std::vector<double> row_upper(b_ub, b_ub + num_row);
    // Left-hand side values (lower bounds)
    std::vector<double> row_lower(num_row, -kHighsInf);

    // Variable bounds
    std::vector<double> col_lower(num_col);
    std::vector<double> col_upper(num_col);

    for (int i = 0; i < num_col; ++i)
    {
        col_lower[i] = bounds[2 * i];
        col_upper[i] = bounds[2 * i + 1];
    }

    // Add columns and rows to HiGHS
    highs.addCols(num_col, col_cost.data(), col_lower.data(), col_upper.data(),
                  0, nullptr, nullptr, nullptr);
    highs.addRows(num_row, row_lower.data(), row_upper.data(),
                  A_value.size(), A_start.data(), A_index.data(), A_value.data());

    // Run HiGHS
    highs.run();

    // Get solution
    HighsSolution solution = highs.getSolution();
    HighsModelStatus model_status = highs.getModelStatus();

    // Prepare the result
    result->x = new double[num_col];
    for (int i = 0; i < num_col; ++i)
    {
        result->x[i] = solution.col_value[i];
    }
    result->fun = highs.getObjectiveValue();
    result->status = static_cast<int>(model_status);
    strncpy(result->message, highs.modelStatusToString(model_status).c_str(), 128);

    return result;
}




std::tuple<std::vector<double>, double, int, std::string>
linprog_highs(const std::vector<double>& c, const std::vector<std::vector<double>>& A_ub,
              const std::vector<double>& b_ub, const std::vector<std::pair<double, double>>& bounds) {

    int num_vars = c.size();
    int num_constraints = b_ub.size();

    // Converti il vettore c in un array di double
    std::vector<double> c_array(c.begin(), c.end());

    // Converti la matrice A_ub in un array piatto di double
    std::vector<double> A_ub_array;
    A_ub_array.reserve(num_constraints * num_vars);
    for (const auto& row : A_ub) {
        A_ub_array.insert(A_ub_array.end(), row.begin(), row.end());
    }

    // Converti il vettore b_ub in un array di double
    std::vector<double> b_ub_array(b_ub.begin(), b_ub.end());

    // Converti i bounds in un array piatto di double
    std::vector<double> bounds_array;
    bounds_array.reserve(2 * num_vars);
    for (const auto& bound : bounds) {
        bounds_array.push_back(bound.first);
        bounds_array.push_back(bound.second);
    }

    // Chiama la funzione C++ linprog_highs con i parametri convertiti
    LinprogResult* result = linprog_highs(c_array.data(), A_ub_array.data(), b_ub_array.data(),
                                          bounds_array.data(), num_vars, num_constraints);

    // Estrai i risultati dalla struttura LinprogResult
    std::vector<double> solution(result->x, result->x + num_vars);
    double fun = result->fun;
    int status = result->status;
    std::string message(result->message);

    // Libera la memoria allocata per il risultato
    free_linprog_result(result);

    // Restituisci i risultati come tuple
    return std::make_tuple(solution, fun, status, message);
}


std::vector<Cell> searchmincells_lp(const QNode& leaf, const std::vector<std::string>& hamstrings)
{
    std::vector<Cell> cells;
    int dims = leaf.getMBR().size();
    std::vector<HalfSpace> leaf_covered = leaf.getCovered();

    auto halfspaces = leaf.getHalfspaces();
    if (halfspaces.empty()) {
        std::vector<std::array<double, 2>> mbr = leaf.getMBR();
        std::vector<double> feasible_coords(mbr.size());
        for (size_t i = 0; i < mbr.size(); ++i) {
            feasible_coords[i] = (mbr[i][0] + mbr[i][1]) / 2.0;
        }
        Point feasible_pnt(feasible_coords);
        std::vector<HalfSpace> empty_halfspaces;
        return {Cell(0, "", leaf_covered, empty_halfspaces, mbr, feasible_pnt)};
    }

    Eigen::VectorXd c = Eigen::VectorXd::Zero(dims + 1);
    c[dims] = -1;

    Eigen::MatrixXd A_ub = Eigen::MatrixXd::Ones(halfspaces.size() + 1, dims + 1);
    A_ub(halfspaces.size(), dims) = 0;

    Eigen::VectorXd b_ub = Eigen::VectorXd::Ones(halfspaces.size() + 1);

    // Configurazione dei bounds
    std::vector<std::pair<double, double>> bounds(dims + 1);
    for (int d = 0; d < dims; ++d) {
        bounds[d] = {leaf.getMBR()[d][0], leaf.getMBR()[d][1]};  // Limiti per le prime dims variabili
    }
    bounds[dims] = {0, std::numeric_limits<double>::infinity()}; // Limite per la variabile slack

    for (const auto& hamstr : hamstrings) {
        for (int b = 0; b < hamstr.size(); ++b) {
            if (hamstr[b] == '0') {
                A_ub.row(b).head(dims) = -halfspaces[b].coeff;
                b_ub[b] = -halfspaces[b].known;
            } else {
                A_ub.row(b).head(dims) = halfspaces[b].coeff;
                b_ub[b] = leaf.getHalfspaces()[b].known;
            }
        }

        std::vector<std::vector<double>> A_ub_std(A_ub.rows(), std::vector<double>(A_ub.cols()));
        for (int i = 0; i < A_ub.rows(); ++i) {
            for (int j = 0; j < A_ub.cols(); ++j) {
                A_ub_std[i][j] = A_ub(i, j);
            }
        }

        std::vector<double> c_std(c.data(), c.data() + c.size());
        std::vector<double> b_ub_std(b_ub.data(), b_ub.data() + b_ub.size());

        auto [solution, fun, status, message] = linprog_highs(c_std, A_ub_std, b_ub_std, bounds);

        if (status == static_cast<int>(HighsModelStatus::kOptimal)) {
            Point feasible_pnt(std::vector<double>(solution.begin(), solution.end() - 1));
            Cell cell(0, hamstr, leaf_covered, leaf.getHalfspaces(), leaf.getMBR(), feasible_pnt);
            cells.push_back(cell);
            break;
        }
    }

    return cells;
}

