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
#include <cstring> // Per memcpy

int numOfSubdivisions = 0;

void freeCell(Cell* cell) {
    if (cell->mask) {
        free(cell->mask);
    }
    if (cell->covered) {
        free(cell->covered);
    }
    if (cell->halfspaces) {
        free(cell->halfspaces);
    }
    if (cell->leaf_mbr) {
        for (int i = 0; i < cell->dims; ++i) {
            free(cell->leaf_mbr[i]);
        }
        free(cell->leaf_mbr);
    }
    // Explicitly call the destructor of feasible_pnt
    cell->feasible_pnt.~Point();

    free(cell);
}

Cell* deepCopyCell(const Cell* original) {
    // Allocate a new Cell
    Cell* copy = (Cell*)malloc(sizeof(Cell));

    // Copy primitive fields
    copy->order = original->order;
    copy->numCovered = original->numCovered;
    copy->numHalfspaces = original->numHalfspaces;
    copy->dims = original->dims;

    // Copy mask
    if (original->mask) {
        copy->mask = (char*)malloc((strlen(original->mask) + 1) * sizeof(char));
        strcpy(copy->mask, original->mask);
    } else {
        copy->mask = nullptr;
    }

    // Copy covered pointers
    if (original->covered && original->numCovered > 0) {
        copy->covered = (HalfSpace**)malloc(original->numCovered * sizeof(HalfSpace*));
        memcpy(copy->covered, original->covered, original->numCovered * sizeof(HalfSpace*));
    } else {
        copy->covered = nullptr;
    }

    // Copy halfspaces pointers
    if (original->halfspaces && original->numHalfspaces > 0) {
        copy->halfspaces = (HalfSpace**)malloc(original->numHalfspaces * sizeof(HalfSpace*));
        memcpy(copy->halfspaces, original->halfspaces, original->numHalfspaces * sizeof(HalfSpace*));
    } else {
        copy->halfspaces = nullptr;
    }

    // Copy leaf_mbr
    if (original->leaf_mbr && original->dims > 0) {
        copy->leaf_mbr = (double**)malloc(original->dims * sizeof(double*));
        for (int i = 0; i < original->dims; ++i) {
            copy->leaf_mbr[i] = (double*)malloc(2 * sizeof(double));
            memcpy(copy->leaf_mbr[i], original->leaf_mbr[i], 2 * sizeof(double));
        }
    } else {
        copy->leaf_mbr = nullptr;
    }

    // Copy feasible_pnt using the copy constructor
    new (&copy->feasible_pnt) Point(original->feasible_pnt);

    return copy;
}



std::pair<int, Cell**> aa_hd(Point** data, int data_size, const Point& p, int& numMinCellsToReturn) {
    numOfSubdivisions = (int)pow(2.0, p.dims - 1);
    QTree qt(p.dims - 1, 10);

    // Ottieni i dominatori
    Point** dominators = nullptr;
    int numDominators = 0;
    getdominators(data, data_size, p, &dominators, numDominators);

    // Ottieni gli incomparabili
    Point** incomp = nullptr;
    int numIncomp = 0;
    getincomparables(data, data_size, p, &incomp, numIncomp);

    // Inizializza una lista per tenere traccia di tutti gli halfspace creati
    HalfSpace** allHalfSpaces = nullptr;
    int numAllHalfSpaces = 0;

    // Lista cumulativa di tutti i record già processati
    Point** all_old_records = nullptr;
    int numAllOldRecords = 0;

    // Definisci la funzione updateqt
    auto updateqt = [&](Point** old_sky, int numOldSky, Point**& new_sky, int& numNewSky, QNode**& leaves, int& numLeaves) {
        // Ottieni il nuovo skyline
        getskyline(incomp, numIncomp, &new_sky, numNewSky);

        // Genera gli halfspaces senza duplicati
        int numNewHalfspaces = 0;
        HalfSpace** new_halfspaces = genhalfspaces(p, new_sky, all_old_records, numNewSky, numAllOldRecords, numNewHalfspaces);

        // Aggiungi gli halfspaces alla lista globale per la gestione della memoria
        for (int i = 0; i < numNewHalfspaces; ++i) {
            numAllHalfSpaces++;
            allHalfSpaces = (HalfSpace**)realloc(allHalfSpaces, numAllHalfSpaces * sizeof(HalfSpace*));
            allHalfSpaces[numAllHalfSpaces - 1] = new_halfspaces[i];
        }

        if (numNewHalfspaces > 0) {
            qt.inserthalfspaces(new_halfspaces, numNewHalfspaces);
            std::cout << "> " << numNewHalfspaces << " halfspace(s) have been inserted" << std::endl;
        }

        free(new_halfspaces); // Libera l'array di puntatori

        // Aggiungi i nuovi punti a `all_old_records`
        numAllOldRecords += numNewSky;
        all_old_records = (Point**)realloc(all_old_records, numAllOldRecords * sizeof(Point*));
        memcpy(all_old_records + numAllOldRecords - numNewSky, new_sky, numNewSky * sizeof(Point*));

        // Ottieni le foglie
        leaves = qt.getleaves(numLeaves);

        for (int i = 0; i < numLeaves; ++i) {
            leaves[i]->setOrder();
        }

        // Ordina le foglie in base all'ordine
        std::sort(leaves, leaves + numLeaves, [](QNode* a, QNode* b) {
        return a->order < b->order;  // Ordina in ordine crescente
    });
    };

    Point** sky = nullptr;
    int numSky = 0;
    QNode** leaves = nullptr;
    int numLeaves = 0;

    updateqt(nullptr, 0, sky, numSky, leaves, numLeaves);

    int minorder_singular = std::numeric_limits<int>::max();
    Cell** mincells_singular = nullptr;
    int n_exp = 0;

    while (true) {
        std::cout << "Cycle number " << n_exp << std::endl;
        int minorder = std::numeric_limits<int>::max();
        Cell** mincells = nullptr;
        int numMinCells = 0;

        for (int i = 0; i < numLeaves; ++i) {
            QNode* leaf = leaves[i];
            int leaf_order = static_cast<int>(leaf->order);
            if (leaf_order > minorder || leaf_order > minorder_singular) {
                break;
            }

            int hamweight = 0;
            while (hamweight <= leaf->numHalfspaces && leaf_order + hamweight <= minorder && leaf_order + hamweight <= minorder_singular) {
                int numHamstrings = 0;
                char** hamstrings = genhammingstrings(static_cast<int>(leaf->numHalfspaces), hamweight, numHamstrings);

                int numCells = 0;
                Cell** cells = searchmincells_lp(*leaf, hamstrings, numHamstrings, numCells);

                // Libera hamstrings
                for (int h = 0; h < numHamstrings; ++h) {
                    free(hamstrings[h]);
                }
                free(hamstrings);

                if (numCells > 0) {
                    for (int c = 0; c < numCells; ++c) {
                        cells[c]->order = leaf_order + hamweight;
                    }

                    if (minorder > leaf_order + hamweight) {
                        minorder = leaf_order + hamweight;
                        // Libera le precedenti mincells
                        if (mincells) {
                            for (int c = 0; c < numMinCells; ++c) {
                                freeCell(mincells[c]);
                            }
                            free(mincells);
                        }

                        // Alloca nuove celle e copia una per una
                        mincells = (Cell**)malloc(numCells * sizeof(Cell*));
                        for (int c = 0; c < numCells; ++c) {
                            mincells[c] = deepCopyCell(cells[c]);  // Copia profonda
                        }
                        numMinCells = numCells;

                    } else {
                        // Aggiungi cells a mincells
                        int totalCells = numMinCells + numCells;
                        Cell** newMincells = (Cell**)malloc(totalCells * sizeof(Cell*));

                        // Copia le celle già presenti in mincells
                        for (int i = 0; i < numMinCells; ++i) {
                            newMincells[i] = deepCopyCell(mincells[i]);  // Copia profonda
                            freeCell(mincells[i]);  // Libera le vecchie celle
                        }
                        free(mincells);

                        // Copia le nuove celle da cells
                        for (int i = 0; i < numCells; ++i) {
                            newMincells[numMinCells + i] = deepCopyCell(cells[i]);  // Copia profonda
                        }

                        mincells = newMincells;
                        numMinCells = totalCells;

                        // Libera cells se non ti serve più
                        for (int c = 0; c < numCells; ++c) {
                            freeCell(cells[c]);  // Libera correttamente
                        }
                        free(cells);
                    }

                    break;
                }
                hamweight++;
            }
        }
        std::cout << "> Expansion " << n_exp << ": Found " << numMinCells << " mincell(s)" << std::endl;

        int new_singulars = 0;
        HalfSpace** to_expand = nullptr;
        int numToExpand = 0;
        for (int c = 0; c < numMinCells; ++c) {
            Cell* cell = mincells[c];
            if (cell->issingular()) {
                minorder_singular = cell->order;
                numMinCellsToReturn++;
                mincells_singular = (Cell**)realloc(mincells_singular, numMinCellsToReturn * sizeof(Cell*));
                mincells_singular[numMinCellsToReturn - 1] = cell;
                /*std::cout << "\n\nMINCELL[" << numMinCellsToReturn - 1 << "]" << std::endl;  // DEBUG
                for (int d = 0; d < mincells_singular[numMinCellsToReturn - 1]->feasible_pnt.dims; ++d) {
                    double coord = mincells_singular[numMinCellsToReturn - 1]->feasible_pnt.coord[d];
                    std::cout << "Coordinate[" << d << "]: " << coord << std::endl;  // DEBUG
                }*/
                new_singulars++;
            } else {
                //std::cout << "\n" << c << "(" << cell->numCovered << ")" << std::endl;
                for (int i = 0; i < cell->numCovered; ++i) {
                    HalfSpace* hs = cell->covered[i];
                    //std::cout << hs->pntID << ",";
                    if (hs->arr == AUGMENTED) {
                        bool found = false;
                        for (int j = 0; j < numToExpand; ++j) {
                            if (to_expand[j] == hs) {
                                found = true;
                                break;
                            }
                        }
                        if (!found) {
                            numToExpand++;
                            to_expand = (HalfSpace**)realloc(to_expand, numToExpand * sizeof(HalfSpace*));
                            to_expand[numToExpand - 1] = hs;
                        }
                    }
                }
            }
        }
        if (new_singulars > 0) {
            std::cout << "> Expansion " << n_exp << ": Found " << new_singulars << " singular mincell(s) with a minorder of " << minorder_singular << std::endl;
        }

        if (numToExpand == 0) {
            // Dealloca tutti gli halfspace creati
            for (int h = 0; h < numAllHalfSpaces; ++h) {
                delete allHalfSpaces[h];
            }
            free(allHalfSpaces);

            // Ritorna il risultato
            return {static_cast<int>(numDominators) + minorder_singular + 1, mincells_singular};
        }

        n_exp++;
        std::cout << "> Expansion " << n_exp << ": " << numToExpand << " halfspace(s) will be expanded" << std::endl;
        for (int h = 0; h < numToExpand; ++h) {
            HalfSpace* hs = to_expand[h];
            hs->arr = SINGULAR;
            // Rimuovi hs->pntID da incomp
            for (int i = 0; i < numIncomp; ++i) {
                if (incomp[i]->id == hs->pntID) {
                    // Rimuovi da incomp
                    for (int j = i; j < numIncomp - 1; ++j) {
                        incomp[j] = incomp[j + 1];
                    }
                    numIncomp--;
                    incomp = (Point**)realloc(incomp, numIncomp * sizeof(Point*));
                    break;
                }
            }
        }
        free(to_expand);

        // Aggiorna qt
        updateqt(sky, numSky, sky, numSky, leaves, numLeaves);
    }
}
