/*************************************/
/********** By Egor Demidov **********/
/******* Textbook Problem 1.4.1 ******/
/*************************************/

#include <iostream>

#include "csv.h"
#include "heat_equation_1d/heat_equation_1d_ghost_points.h"
#include "heat_equation_1d/heat_equation_1d_one_sided.h"

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

int main() {
    constexpr long M = 10;
    constexpr double dt = 0.004;
    constexpr double nu = 1.0;
    constexpr double t_points[] = {0.06, 0.1, 0.9};
    long t_point_n_values[std::size(t_points)];
    std::ranges::transform(t_points, std::begin(t_point_n_values), [](double t) {
        return static_cast<long>(t / dt);
    });

    auto exact_solution_fun = [nu](double x, double t) {
        return cos(M_PI * x / 2.0) * exp(-nu * M_PI * M_PI / 4.0 * t);
    };

    HeatEquation1D_GhostPoints<BC_L, BC_R, IC> heat_equation_ghost_pts(M, nu, dt);
    HeatEquation1D_OneSided<BC_L, BC_R, IC> heat_equation_one_sided(M, nu, dt);

    std::vector<double> ghost_pt_solutions[std::size(t_points) + 1],
                        one_sided_solutions[std::size(t_points) + 1],
                        exact_solutions[std::size(t_points) + 1];

    auto target_time_step_itr = std::begin(t_point_n_values);
    auto ghost_pt_solutions_itr = std::begin(ghost_pt_solutions);
    auto one_sided_solutions_itr = std::begin(one_sided_solutions);
    auto exact_solutions_itr = std::begin(exact_solutions);

    *ghost_pt_solutions_itr = std::vector<double>(heat_equation_ghost_pts.get_u0().begin(), heat_equation_ghost_pts.get_u0().end());
    *one_sided_solutions_itr = std::vector<double>(heat_equation_one_sided.get_u0().begin(), heat_equation_one_sided.get_u0().end());
    ghost_pt_solutions_itr ++;
    one_sided_solutions_itr ++;

    for (long n = 1; n <= t_point_n_values[std::size(t_point_n_values) - 1]; n ++) {
        auto sol_one_sided = heat_equation_one_sided.do_step();
        auto sol_ghost = heat_equation_ghost_pts.do_step();

        if (n == *target_time_step_itr) {
            *ghost_pt_solutions_itr = std::vector<double>(sol_ghost.begin(), sol_ghost.end());
            *one_sided_solutions_itr = std::vector<double>(sol_one_sided.begin(), sol_one_sided.end());

            target_time_step_itr ++;
            ghost_pt_solutions_itr ++;
            one_sided_solutions_itr ++;
        }
    }

    // Populate the exact solution buffer
    std::vector<double> t_points_exact_sol(std::begin(t_points), std::end(t_points));
    t_points_exact_sol.insert(t_points_exact_sol.begin(), 0.0);
    for (double t : t_points_exact_sol) {
        exact_solutions_itr->resize(M + 1);
        for (long m = 0; m < M + 1; m ++) {
            double xm = 1.0 / static_cast<double>(M) * static_cast<double>(m);
            (*exact_solutions_itr)[m] = exact_solution_fun(xm, t);
        }
        exact_solutions_itr ++;
    }

    CSVWriter<5> ghost_pt_writer("sol_ghost_point.csv", {"x", "t=0", "t=0.06", "t=0.1", "t=0.9"});
    CSVWriter<5> one_sided_writer("sol_one_sided.csv", {"x", "t=0", "t=0.06", "t=0.1", "t=0.9"});
    CSVWriter<5> exact_writer("sol_exact.csv", {"x", "t=0", "t=0.06", "t=0.1", "t=0.9"});

    for (long m = 0; m < M + 1; m ++) {
        double xm = 1.0 / static_cast<double>(M) * static_cast<double>(m);

        std::array<double, std::size(t_points) + 2> ghost_pt_line,
                                                    one_sided_line,
                                                    exact_line;

        ghost_pt_line[0] = xm;
        one_sided_line[0] = xm;
        exact_line[0] = xm;
        for (size_t i = 0; i < std::size(t_points) + 1; i ++) {
            ghost_pt_line[i + 1] = ghost_pt_solutions[i][m];
            one_sided_line[i + 1] = one_sided_solutions[i][m];
            exact_line[i + 1] = exact_solutions[i][m];
        }

        ghost_pt_writer.append_line(ghost_pt_line);
        one_sided_writer.append_line(one_sided_line);
        exact_writer.append_line(exact_line);
    }

    return 0;
}
