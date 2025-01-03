//
// Created by leona on 06/08/2024.
//

#include "maxrank.h"
#include "query.h"
#include "halfspace.h"
#include "cell.h"
#include "qtree.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string.h>

int numOfSubdivisions = 0;

// Definisci e inizializza le variabili globali
HalfSpaceCache* halfspaceCache = nullptr;
std::unordered_map<Point, long, PointHash> pointToHalfSpaceCache;

std::vector<std::string> readCombinations(const int& dims) {
    std::vector<std::string> comb;
    FILE *fp;
    char *token;
    char m_separator[] = " \n\t";
    char buf[512];
    std::string FileName[] = {"../bin/Comb2D.txt", "../bin/Comb3D.txt", "../bin/Comb4D.txt", "../bin/Comb5D.txt", "../bin/Comb6D.txt", "../bin/Comb7D.txt",
                         "../bin/Comb8D.txt", "../bin/Comb9D.txt"};
    std::string strComb;
    fp = fopen(FileName[dims - 2].c_str(), "r");
    if (fp == NULL) {
        std::cout << "error in fileopen!" << std::endl;
        exit(0);
    }
    fgets(buf, 512, fp);
    if (atoi(buf) != dims) {
        std::cout << "Error! Dimensions are not equal!" << std::endl;
        exit(0);
    }
    while (fgets(buf, 512, fp) != NULL) {
        token = strtok(buf, m_separator);
        //while (token != NULL){
        //	token = strtok(NULL,m_separator);
        //}
        std::string strComb = token;
        comb.push_back(strComb);
    }
    fclose(fp);
    return comb;
}
int dims;
float queryPlane[10];
std::vector<std::string> Comb;

bool MbrIsValid(const std::vector<std::array<double, 2>> &mbr) {
    int numAbove = 0;
    int numBelow = 0;
    long int numOfVertices = Comb.size();

    // Usa un array statico per coord invece di ricostruirlo dinamicamente ogni volta
    std::vector<double> coord(dims, 0.0);

    for (const auto &combination : Comb) {
        float sum = 0.0; // Calcola sum durante la costruzione di coord
        for (int j = 0; j < dims; ++j) {
            coord[j] = (combination[j] == '0') ? mbr[j][0] : mbr[j][1];
            sum += coord[j];
        }

        // Confronta sum con hs[dims]
        if (sum > queryPlane[dims]) {
            ++numAbove;
            if (numAbove > 0 && numBelow > 0) break; // Early exit: l'MBR è intersecato
        } else if (sum < queryPlane[dims]) {
            ++numBelow;
            if (numAbove > 0 && numBelow > 0) break; // Early exit: l'MBR è intersecato
        }
    }

    // Restituisci il risultato in base ai conteggi
    if (numAbove == numOfVertices) return false;
    if (numBelow == numOfVertices) return true;
    return false;
}

std::pair<int, std::vector<Cell>> aa_hd(const std::vector<Point>& data, const Point& p) {
    dims = static_cast<int>(p.dims - 1);
    numOfSubdivisions = (int) pow(2.0, dims);
    Comb = readCombinations(dims);
    for (int i = 0; i < dims + 1; i++) queryPlane[i] = 1;

    QTree qt(dims, 10);
    std::vector<Point> dominators = getdominators(data, p);
    std::vector<Point> incomp = getincomparables(data, p);

    // Inizializzo la cache per gli halfspaces
    initializeCache(data.size());

    auto updateqt = [&](const std::vector<Point>& old_sky) {
        std::vector<Point> new_sky = getskyline(incomp);
        std::vector<long> new_halfspaces = genhalfspaces(p, new_sky);
        std::vector<long> unique_new_halfspaces;
        for (const auto& hs : new_halfspaces) {
            bool found = false;
            for (const auto& os : old_sky) {
                if (os.id == hs) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                unique_new_halfspaces.push_back(hs);
            }
        }

        new_halfspaces = std::move(unique_new_halfspaces);

        if (!new_halfspaces.empty()) {
            //qt.inserthalfspaces(new_halfspaces);
            qt.inserthalfspacesMacroSplit(new_halfspaces);
            std::cout << "> " << new_halfspaces.size() << " halfspace(s) have been inserted" << std::endl;
        }

        auto new_leaves = qt.getAllLeaves();//qt.getLeaves();
        qt.updateAllOrders();

        std::sort(new_leaves.begin(), new_leaves.end(), [](QNode* a, QNode* b) { return a->getOrder() < b->getOrder(); });

        return std::make_pair(new_sky, new_leaves);
    };


    auto [sky, leaves] = updateqt({});

    int minorder_singular = std::numeric_limits<int>::max();
    std::vector<Cell> mincells_singular;
    int n_exp = 0;

    while (true) {
        std::cout << "Cycle number " << n_exp << std::endl;
        int minorder = std::numeric_limits<int>::max();
        std::vector<Cell> mincells;

        for (auto leaf : leaves) {
            int leaf_order = static_cast<int>(leaf->getOrder());
            if (leaf_order > minorder || leaf_order > minorder_singular) {
                break;
            }
            //prune away leaf nodes that lie about hyperplane q_1+q2+...+q_d < 1;
            if (!MbrIsValid(leaf->getMBR())) {
                //std::cout << "Leaf node " << leaf->getNodeID() << " is pruned!" << std::endl;
                continue;
            }

            int hamweight = 0;
            while (hamweight <= leaf->getHalfspaces().size() && leaf_order + hamweight <= minorder && leaf_order + hamweight <= minorder_singular) {
                std::vector<std::string> hamstrings = genhammingstrings(static_cast<int>(leaf->getHalfspaces().size()), hamweight);
                std::vector<Cell> cells = searchmincells_lp(*leaf, hamstrings);

                if (!cells.empty()) {
                    for (auto& cell : cells) {
                        cell.order = leaf_order + hamweight;
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
        for (auto& cell : mincells) {
            if (cell.issingular()) {
                minorder_singular = cell.order;
                mincells_singular.push_back(cell);
                new_singulars++;
            } else {
                for (auto k : cell.covered) {
                    auto hs = halfspaceCache->get(k);
                    if (hs->arr == Arrangement::AUGMENTED && std::find(to_expand.begin(), to_expand.end(), hs) == to_expand.end()) {
                        to_expand.push_back(hs);
                    }
                }
            }
        }
        if (new_singulars > 0) {
            std::cout << "> Expansion " << n_exp << ": Found " << new_singulars << " singular mincell(s) with a minorder of " << minorder_singular << std::endl;
        }

        if (to_expand.empty()) {
            return {static_cast<int>(dominators.size()) + minorder_singular + 1, mincells_singular};
        }

        n_exp++;
        std::cout << "> Expansion " << n_exp << ": " << to_expand.size() << " halfspace(s) will be expanded" << std::endl;
        for (const auto& hs : to_expand) {
            hs->arr = Arrangement::SINGULAR;
            auto it = std::find_if(incomp.begin(), incomp.end(), [&](const Point& pt) { return hs->pntID == pt.id; });
            if (it != incomp.end()) {
                incomp.erase(it);
            }

        }
        std::tie(sky, leaves) = updateqt(sky);
    }
}
