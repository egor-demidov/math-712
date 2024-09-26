/*************************************/
/********** By Egor Demidov **********/
/******* Textbook Problem 1.4.1 ******/
/*************************************/

#include <iostream>

#include "csv.h"
#include "heat_equation_1d/heat_equation_1d_ghost_points.h"

// Define the left boundary condition
struct BC_L {
    static constexpr BCType bc_type = NEUMANN;
    static constexpr double evaluate(double t [[maybe_unused]]) {
        return 0.0;
    }
};

// Define the right boundary condition
struct BC_R {
    static constexpr BCType bc_type = DIRICHLET;
    static constexpr double evaluate(double t [[maybe_unused]]) {
        return 0.0;
    }
};

// Define the initial condition
struct IC {
    static double evaluate(double x) {
        return cos(M_PI * x / 2.0);
    }
};

double l2_norm(std::span<const double> data_span) {
    double accumulator = 0.0;
    for (double val : data_span)
        accumulator += val * val;
    return sqrt(accumulator);
}

// Discretization parameters:
// M - number of space intervals (number of points: M+1)
// N - number of time steps
// t_final - final time
// dt - integration time step

int main() {
    constexpr long M = 10;
    constexpr long N = 10;
    constexpr double t_final = 0.1;

    HeatEquation1D<BC_L, BC_R, IC> heat_equation(M, 1.0,  t_final / static_cast<double>(N));

    for (long n = 1; n <= N; n ++) {
        auto sol = heat_equation.do_step();

        std::cout << l2_norm(sol) << "\t\t";
        for (auto u : sol)
            std::cout << u << " ";
        std::cout << "\n";
    }

    return 0;
}
