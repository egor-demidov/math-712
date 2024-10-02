//
// Created by egor on 9/26/24.
//

#include "heat_equation_1d/heat_equation_1d_ghost_points.h"
#include "heat_equation_1d/heat_equation_1d_one_sided.h"
#include "csv.h"

/*
0.1	0.00224249	2.32898e-05
0.01	0.000247977	2.82819e-07
0.001	2.50792e-05	2.20609e-09
0.0001	2.51018e-06	3.71175e-11
*/

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

std::pair<double, double> get_l2_norms_of_error(const long M) {
    constexpr double dt_over_h2 = 0.004 / 0.01;
    constexpr double nu = 1.0;
    constexpr double t_tot = 2.0;
    const double dt = dt_over_h2 / static_cast<double>(M) / static_cast<double>(M);
    const long N = static_cast<long>(t_tot / dt);
    const double h = 1.0 / static_cast<double>(M);

    auto exact_solution_fun = [](double x, double t) {
        return cos(M_PI * x / 2.0) * exp(-nu * M_PI * M_PI / 4.0 * t);
    };

    HeatEquation1D_GhostPoints<BC_L, BC_R, IC> heat_equation_ghost_pts(M, nu, dt);
    HeatEquation1D_OneSided<BC_L, BC_R, IC> heat_equation_one_sided(M, nu, dt);

    std::vector<double> ghost_pt_solution,
            one_sided_solution,
            exact_solution;

    for (long n = 1; n < N - 1; n ++) {
        heat_equation_ghost_pts.do_step();
    }
    auto ghost_pt_solution_span = heat_equation_ghost_pts.do_step();
    ghost_pt_solution = std::vector(ghost_pt_solution_span.begin(), ghost_pt_solution_span.end());

    for (long n = 1; n < N - 1; n ++) {
        heat_equation_one_sided.do_step();
    }
    auto one_sided_solution_span = heat_equation_one_sided.do_step();
    one_sided_solution = std::vector(one_sided_solution_span.begin(), one_sided_solution_span.end());

    exact_solution.resize(M + 1);
    for (long m = 0; m < M + 1; m ++) {
        double xm = 1.0 / static_cast<double>(M) * static_cast<double>(m);
        exact_solution[m] = exact_solution_fun(xm, t_tot);
    }

    double l2_norm_error_ghost_pts = 0.0,
            l2_norm_error_one_sided = 0.0;

    for (long m = 0; m < M + 1; m ++) {
        l2_norm_error_ghost_pts += pow(ghost_pt_solution[m] - exact_solution[m], 2.0) * h;
        l2_norm_error_one_sided += pow(one_sided_solution[m] - exact_solution[m], 2.0) * h;
    }

    return std::make_pair(sqrt(l2_norm_error_one_sided), sqrt(l2_norm_error_ghost_pts));
}

int main() {


    constexpr long M_values[] = {10, 100, 1'000, 10'000, 100'000};

    CSVWriter<3> one_sided_writer("sol_one_sided_part2.csv", {"h", "1st", "2nd"});

    for (long M : M_values) {
        auto [l2_1st_order, l2_2nd_order] = get_l2_norms_of_error(M);
        double h = 1.0 / static_cast<double>(M);
        std::cout << h << "\t" << l2_1st_order << "\t" << l2_2nd_order << "\n";
        one_sided_writer.append_line({h, l2_1st_order, l2_2nd_order});
    }

    return 0;
}
