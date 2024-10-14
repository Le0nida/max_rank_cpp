//
// Created by leona on 06/08/2024.
//

#include "maxrank.h"
#include <algorithm>
#include "query.h"
#include "halfspace.h"
#include "cell.h"
#include "qtree.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <unordered_set>

int numOfSubdivisions = 0;
std::unordered_set<long int> HalfSpaces;

std::pair<int, std::vector<std::shared_ptr<Cell>>> aa_hd(const std::vector<std::shared_ptr<Point>>& data, const Point& p) {
    numOfSubdivisions = static_cast<int>(std::pow(2.0, p.dims - 1));
    QTree qt(p.dims - 1, 10);

    // Ottieni i dominatori
    std::vector<std::shared_ptr<Point>> dominators;
    getdominators(data, p, dominators);

    // Ottieni gli incomparabili
    std::vector<std::shared_ptr<Point>> incomp;
    getincomparables(data, p, incomp);

    std::vector<std::shared_ptr<HalfSpace>> halfspacesToInsert;

    // Funzione lambda per aggiornare il QTree
    auto updateqt = [&](std::vector<std::shared_ptr<Point>>& new_sky, std::vector<QNode*>& leaves) {
        getskyline(incomp, new_sky);

        halfspacesToInsert.clear();
        halfspacesToInsert = genhalfspaces(p, new_sky);

        if (!halfspacesToInsert.empty()) {
            qt.inserthalfspaces(halfspacesToInsert);
            std::cout << "> " << halfspacesToInsert.size() << " halfspace(s) have been inserted" << std::endl;
        }

        leaves = qt.getleaves();

        for (const auto& leaf : leaves) {
            leaf->setOrder();
        }

        std::sort(leaves.begin(), leaves.end(), [](const QNode* a, const QNode* b) {
            return a->order < b->order;
        });
    };

    std::vector<std::shared_ptr<Point>> sky;
    std::vector<QNode*> leaves;
    updateqt(sky, leaves);

    int minorder_singular = std::numeric_limits<int>::max();
    std::vector<std::shared_ptr<Cell>> mincells_singular;
    int n_exp = 0;

    while (true) {
        std::cout << "Cycle number " << n_exp << std::endl;
        int minorder = std::numeric_limits<int>::max();
        std::vector<std::shared_ptr<Cell>> mincells;

        for (auto& leaf : leaves) {
            int leaf_order = static_cast<int>(leaf->order);
            if (leaf_order > minorder || leaf_order > minorder_singular) {
                break;
            }

            int hamweight = 0;
            while (hamweight <= leaf->halfspaces.size() && leaf_order + hamweight <= minorder && leaf_order + hamweight <= minorder_singular) {
                int numHamstrings = 0;
                std::unique_ptr<char*[]> hamstrings = std::unique_ptr<char*[]>(genhammingstrings(static_cast<int>(leaf->halfspaces.size()), hamweight, numHamstrings));

                std::vector<std::shared_ptr<Cell>> cells = searchmincells_lp(*leaf, hamstrings.get(), numHamstrings);
                if (!cells.empty()) {
                    for (auto& cell : cells) {
                        cell->order = leaf_order + hamweight;
                    }

                    if (minorder > leaf_order + hamweight) {
                        minorder = leaf_order + hamweight;
                        mincells = std::move(cells);
                    } else {
                        mincells.insert(mincells.end(), cells.begin(), cells.end());
                    }
                    break;
                }
                hamweight++;
            }
        }

        std::cout << "> Expansion " << n_exp << ": Found " << mincells.size() << " mincell(s)" << std::endl;

        int new_singulars = 0;
        std::vector<std::shared_ptr<HalfSpace>> to_expand;
        int i = 0;
        for (auto& cell : mincells) {
            if (cell->issingular()) {
                minorder_singular = cell->order;
                mincells_singular.push_back(cell);
                new_singulars++;
            } else {
                for (auto& hs : cell->covered) {
                    if (hs->arr == AUGMENTED && std::find(to_expand.begin(), to_expand.end(), hs) == to_expand.end()) {
                        to_expand.push_back(hs);
                    }
                }
            }
        }

        if (new_singulars > 0) {
            std::cout << "> Expansion " << n_exp << ": Found " << new_singulars << " singular mincell(s) with a minorder of " << minorder_singular << std::endl;
        }

        if (to_expand.empty()) {
            return {static_cast<int>(dominators.size()) + minorder_singular + 1, std::move(mincells_singular)};
        }

        n_exp++;
        std::cout << "> Expansion " << n_exp << ": " << to_expand.size() << " halfspace(s) will be expanded" << std::endl;

        for (const auto& hs : to_expand) {
            hs->arr = SINGULAR;
            auto it = std::remove_if(incomp.begin(), incomp.end(), [&](const std::shared_ptr<Point>& point) {
                return point->id == hs->pntID;
            });
            incomp.erase(it, incomp.end());
        }

        updateqt(sky, leaves);
    }
}
