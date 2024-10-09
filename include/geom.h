//
// Created by leona on 01/08/2024.
//

#ifndef GEOM_H
#define GEOM_H

#include <cstdlib> // Per malloc e free

class Point {
public:
    // Costruttore
    Point(double* coord, int dims, int id = -1);
    // Distruttore
    ~Point();

    Point(const Point& other);
    Point& operator=(const Point& other);

    // Proprietà pubbliche
    int id;
    double* coord;
    int dims;

    // Operatore di uguaglianza per confrontare due punti
    bool operator==(const Point& other) const;
};

// Funzione per generare le maschere
void genmasks(int dims, double***& pts_mask, int& num_pts, double***& nds_mask, int& num_nds);

#endif // GEOM_H

