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

std::vector<Cell> searchmincells_lp(const QNode& leaf, const std::vector<std::string>& hamstrings)
{
    std::vector<Cell> cells;
    int dims = leaf.getMBR().size();
    std::vector<HalfSpace> leaf_covered = leaf.getCovered();

    // If there are no halfspaces, then the whole leaf is the mincell
    auto halfspaces = leaf.getHalfspaces();
    if (halfspaces.empty())
    {
        std::vector<std::array<double, 2>> mbr = leaf.getMBR();
        std::vector<double> feasible_coords(mbr.size());
        for (size_t i = 0; i < mbr.size(); ++i)
        {
            feasible_coords[i] = (mbr[i][0] + mbr[i][1]) / 2.0; // Scegli un punto centrale nel bounding box
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

    std::vector<std::pair<double, double>> bounds(dims);
    for (int d = 0; d < dims; ++d)
    {
        bounds[d] = {leaf.getMBR()[d][0], leaf.getMBR()[d][1]};
    }
    bounds.push_back({0, std::numeric_limits<double>::infinity()});

    for (const auto& hamstr : hamstrings)
    {
        for (int b = 0; b < hamstr.size(); ++b)
        {
            if (hamstr[b] == '0')
            {
                A_ub.row(b).head(dims) = -halfspaces[b].coeff;
                b_ub[b] = -halfspaces[b].known;
            }
            else
            {
                A_ub.row(b).head(dims) = halfspaces[b].coeff;
                b_ub[b] = halfspaces[b].known;
            }
        }

        // Convert to appropriate types for the linear programming solver
        std::vector<double> c_std(c.data(), c.data() + c.size());
        std::vector<std::vector<double>> A_ub_std;
        for (int i = 0; i < A_ub.rows(); ++i)
        {
            A_ub_std.push_back(std::vector<double>(A_ub.row(i).data(), A_ub.row(i).data() + A_ub.cols()));
        }
        std::vector<double> b_ub_std(b_ub.data(), b_ub.data() + b_ub.size());


        Highs highs;
        HighsLp model;

        // Configurazione di base del modello LP
        model.num_col_ = dims + 1;
        model.num_row_ = A_ub_std.size();
        model.sense_ = ObjSense::kMaximize;

        // Definizione dei costi della funzione obiettivo
        model.col_cost_ = c_std;

        // Definizione dei limiti sulle variabili
        model.col_lower_ = std::vector<double>(dims + 1, 0);
        model.col_upper_ = std::vector<double>(dims + 1, std::numeric_limits<double>::infinity());

        // Definizione dei limiti sulle righe
        model.row_lower_ = std::vector<double>(A_ub_std.size(), -std::numeric_limits<double>::infinity());
        model.row_upper_ = b_ub_std;

        // Configurazione della matrice sparse A
        HighsSparseMatrix matrix;
        matrix.num_col_ = dims + 1;
        matrix.num_row_ = A_ub_std.size();

        std::vector<HighsInt> start(A_ub_std.size() + 1, 0);
        std::vector<HighsInt> index;
        std::vector<double> value;

        int index_counter = 0;
        for (int i = 0; i < A_ub_std.size(); ++i)
        {
            start[i] = index_counter;
            for (int j = 0; j < dims + 1; ++j)
            {
                if (A_ub_std[i][j] != 0.0)
                {
                    // Aggiungi solo i valori non zero
                    index.push_back(j);
                    value.push_back(A_ub_std[i][j]);
                    ++index_counter;
                }
            }
        }
        start[A_ub_std.size()] = index_counter;

        matrix.start_ = std::move(start);
        matrix.index_ = std::move(index);
        matrix.value_ = std::move(value);

        // Assegna la matrice al modello
        model.a_matrix_ = matrix;

        // Carica il modello nel solver
        HighsStatus status = highs.passModel(model);

        // Esegui l'ottimizzazione
        if (status == HighsStatus::kOk)
        {
            status = highs.run();

            if (status == HighsStatus::kOk)
            {
                const auto& highs_solution = highs.getSolution();
                std::vector<double> solution(highs_solution.col_value.begin(), highs_solution.col_value.begin() + dims);

                Point feasible_pnt(solution);
                Cell cell(0, hamstr, leaf_covered, leaf.getHalfspaces(), leaf.getMBR(), feasible_pnt);
                cells.push_back(cell);
            }
            else
            {
                std::cerr << "Highs solver failed to run." << std::endl;
            }
        }
        else
        {
            std::cerr << "Failed to pass model to Highs solver." << std::endl;
        }
    }

    return cells;
}
