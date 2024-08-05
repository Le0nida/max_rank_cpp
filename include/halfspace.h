//
// Created by leona on 01/08/2024.
//

#ifndef HALFSPACE_H
#define HALFSPACE_H

#include <vector>

class Halfspace {
public:
    std::vector<double> coeff;
    double known;
};

enum Position {
    IN = 1,
    OUT = -1,
    on = 0
};

#endif //HALFSPACE_H
