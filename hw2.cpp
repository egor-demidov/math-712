/*************************************/
/********** By Egor Demidov **********/
/******* Textbook Problem 1.4.1 ******/
/*************************************/

#include <iostream>
#include <numeric>

#include <Eigen/Eigen>

#include "csv.h"

// enum BCType {
//     DIRICHLET, NEUMANN
// };
//
// using GridPoints = Eigen::Vector<double, Eigen::Dynamic>;
//
// BC_L, BC_R are functors that compute boundary conditions as a function of time
// template <typename BC_L, typename BC_R, typename IC>
// class HeatEquation1D {
// public:
//     HeatEquation1D(long init_m, double init_nu, GridPoints init_u0, double init_h, double init_dt)
//         : m{init_m}
//         , nu{init_nu}
//         , u0{std::move(init_u0)}
//         , u_size{u0.size()}
//         , h{init_h}
//         , dt{init_dt}
//     {
//         if constexpr (BC_L::type == NEUMANN) {
//             u_size += 1;
//         }
//         if constexpr (BC_R::type == NEUMANN) {
//             u_size += 1;
//         }
//     }
//
// private:
//     long m;
//     double nu;
//     GridPoints u0;
//     GridPoints u;
//     long u_size;
//     double h, dt;
// };
//
// struct BC_L {
//     static constexpr BCType bc_type = DIRICHLET;
//     static constexpr double evaluate(double t [[maybe_unused]]) {
//         return 0.0;
//     }
// };
//
// struct BC_R {
//     static constexpr BCType bc_type = NEUMANN;
//     static constexpr double evaluate(double t [[maybe_unused]]) {
//         return 0.0;
//     }
// };
//
// struct IC {
//     static double evaluate(double x) {
//         return cos(M_PI * x / 2.0);
//     }
// };

// Discretization parameters:
// M - number of space intervals (number of points: M+1)
// N - number of time steps
// t_final - final time
// dt - integration time step

int main() {
    // HeatEquation1D<BC_L, BC_R, IC> heat_equation{};

    // Diffusion coefficient
    constexpr double NU = 1.0;

    // Neumann BC on the left boundary
    auto v_x_l = [](double t) -> double {
        return 0;
    };

    // Dirichlet BC on the right boundary
    auto v_r = [](double t) -> double {
        return 0;
    };

    // Initial condition
    auto v0 = [](double x) -> double {
        return cos(M_PI * x / 2.0);
    };

    constexpr unsigned long M = 10;
    constexpr unsigned long N = 100;
    constexpr double T_FINAL = 0.01;
    constexpr double H = 1.0 / static_cast<double>(M);
    constexpr double DT = T_FINAL / static_cast<double>(N);

    // M+1 points plus an extra ghost point for the Neumann BC on the left
    constexpr unsigned long N_X_POINTS = M + 2;

    Eigen::Vector<double, N_X_POINTS> v0_grid, v_grid;
    for (unsigned long m = 0; m < M + 1; m ++) {
        double x_m = H * static_cast<double>(m);
        v0_grid(m+1) = v0(x_m);
    }
    v0_grid(0) = v0_grid(2) - v_x_l(0) * 2.0 * H;  // Value of the ghost point is prescribed by the Neumann B.C.

    // Do the time steps
    v_grid = v0_grid;
    for (unsigned long n = 1; n < N; n ++) {
        Eigen::Vector<double, N_X_POINTS> v_new_grid;

        double t_new = static_cast<double>(n) * DT;

        // Apply the Dirichlet BC
        v_new_grid(M+1) = v_r(t_new);

        // Skipping point M because of Dirichlet BC
        for (unsigned long m = 1; m < M; m ++) {
            v_new_grid(m+1) = v_grid(m+1) + NU * DT / H / H * (v_grid(m+2) - 2.0 * v_grid(m+1) + v_grid(m));
        }
        v_new_grid(0) = v_new_grid(2) - v_x_l(t_new) * 2.0 * H;  // Value of the ghost point is prescribed by the Neumann B.C.
        v_grid = std::move(v_new_grid);
    }

    CSVWriter<3> csv_writer("sol_ghost_point.csv", {"x", "v0", "v"});
    for (unsigned long m = 1; m < M + 2; m ++) {
        csv_writer.append_line({static_cast<double>(m-1) * H, v0_grid(m), v_grid(m)});
    }

    std::cout << "Ghost point at step 0: " << v0_grid(0) << std::endl;

    return 0;
}
