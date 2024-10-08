#include "query.h"
#include <cstdlib> // Per malloc, free

void getdominators(Point** data, int data_size, const Point& p, Point*** dominators, int& num_dominators) {
    num_dominators = 0;
    *dominators = nullptr;

    for (int idx = 0; idx < data_size; ++idx) {
        Point* r = data[idx];
        bool less_equal = true;
        bool less = false;

        for (int i = 0; i < p.dims; ++i) {
            if (r->coord[i] > p.coord[i]) {
                less_equal = false;
                break;
            }
            if (r->coord[i] < p.coord[i]) {
                less = true;
            }
        }

        if (less_equal && less) {
            num_dominators++;
            *dominators = (Point**)realloc(*dominators, num_dominators * sizeof(Point*));
            (*dominators)[num_dominators - 1] = r;
        }
    }
}

void getdominees(Point** data, int data_size, const Point& p, Point*** dominees, int& num_dominees) {
    num_dominees = 0;
    *dominees = nullptr;

    for (int idx = 0; idx < data_size; ++idx) {
        Point* r = data[idx];
        bool greater_equal = true;
        bool greater = false;

        for (int i = 0; i < p.dims; ++i) {
            if (r->coord[i] < p.coord[i]) {
                greater_equal = false;
                break;
            }
            if (r->coord[i] > p.coord[i]) {
                greater = true;
            }
        }

        if (greater_equal && greater) {
            num_dominees++;
            *dominees = (Point**)realloc(*dominees, num_dominees * sizeof(Point*));
            (*dominees)[num_dominees - 1] = r;
        }
    }
}

void getincomparables(Point** data, int data_size, const Point& p, Point*** incomparables, int& num_incomparables) {
    num_incomparables = 0;
    *incomparables = nullptr;

    for (int idx = 0; idx < data_size; ++idx) {
        Point* r = data[idx];
        bool less = false;
        bool greater = false;

        for (int i = 0; i < p.dims; ++i) {
            if (r->coord[i] < p.coord[i]) {
                less = true;
            }
            if (r->coord[i] > p.coord[i]) {
                greater = true;
            }
        }

        if (less && greater) {
            num_incomparables++;
            *incomparables = (Point**)realloc(*incomparables, num_incomparables * sizeof(Point*));
            (*incomparables)[num_incomparables - 1] = r;
        }
    }
}

void getskyline(Point** data, int data_size, Point*** skyline, int& num_skyline) {
    num_skyline = 0;
    *skyline = nullptr;

    auto dominates = [](const Point* p, const Point* r) {
        bool less_equal = true;
        bool less = false;

        for (int i = 0; i < p->dims; ++i) {
            if (p->coord[i] > r->coord[i]) {
                less_equal = false;
                break;
            }
            if (p->coord[i] < r->coord[i]) {
                less = true;
            }
        }

        return less_equal && less;
    };

    for (int idx = 0; idx < data_size; ++idx) {
        Point* pnt = data[idx];
        bool dominated = false;

        for (int w_idx = 0; w_idx < num_skyline; ++w_idx) {
            if (dominates((*skyline)[w_idx], pnt)) {
                dominated = true;
                break;
            }
        }

        if (!dominated) {
            // Rimuovi i punti nel window che sono dominati da pnt
            int new_num_skyline = 0;
            Point** new_skyline = nullptr;
            for (int w_idx = 0; w_idx < num_skyline; ++w_idx) {
                if (!dominates(pnt, (*skyline)[w_idx])) {
                    new_num_skyline++;
                    new_skyline = (Point**)realloc(new_skyline, new_num_skyline * sizeof(Point*));
                    new_skyline[new_num_skyline - 1] = (*skyline)[w_idx];
                }
            }
            free(*skyline);
            *skyline = new_skyline;
            num_skyline = new_num_skyline;

            // Aggiungi pnt al window
            num_skyline++;
            *skyline = (Point**)realloc(*skyline, num_skyline * sizeof(Point*));
            (*skyline)[num_skyline - 1] = pnt;
        }
    }
}
