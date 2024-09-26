//
// Created by egor on 9/26/24.
//

#ifndef COMMON_DEFS_H
#define COMMON_DEFS_H

#include <vector>
#include <algorithm>
#include <ranges>
#include <cmath>
#include <numeric>

enum BCType {
    DIRICHLET, NEUMANN
};

using GridPoints = std::vector<double>;

#endif //COMMON_DEFS_H
