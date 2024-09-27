//
// Created by egor on 9/26/24.
//

#ifndef HEAT_EQUATION_1D_ONE_SIDED_H
#define HEAT_EQUATION_1D_ONE_SIDED_H

#include "../common_defs.h"

// In this version, Neumann BCs are implemented using one-sided differences

// BC_L, BC_R are functors that compute boundary conditions as a function of time
template <typename BC_L, typename BC_R, typename IC>
class HeatEquation1D_OneSided {
public:
    HeatEquation1D_OneSided(long init_M, double init_nu, double init_dt)
        : M{init_M}
        , nu{init_nu}
        , t{0.0}
        , h{1.0 / M}
        , dt{init_dt}
    {
        // Figure out the buffer size with ghost points
        long u_size = M + 1;

        // Allocate and initialize buffers
        u.resize(u_size);
        u0.resize(u_size);
        u_new.resize(u_size);
        std::ranges::fill(u, 0.0);
        std::ranges::fill(u_new, 0.0);
        std::ranges::fill(u0, 0.0);

        // Initialize the initial condition grid function
        for (long m = 0; m < static_cast<long>(u0.size()); m ++) {
            double xm = static_cast<double>(m) * h;
            u0[m] = IC::evaluate(xm);
        }
        apply_bcs(u0);

        // Copy the initial condition to the current step grid function buffer
        u = u0;
    }

    [[nodiscard]]
    std::span<const double> get_u0() const {
        return std::span<const double>(u);
    }

    std::span<const double> do_step() {
        #pragma omp parallel for
        for (long m = 1; m < static_cast<long>(u.size()) - 1; m ++) {
            // To manipulate array data, use index i which is equivalent to theoretical index m
            u_new[m] = u[m] + nu * dt / h / h * (u[m+1] - 2.0 * u[m] + u[m-1]);
        }

        t += dt;

        apply_bcs(u_new);

        auto temp_vec = std::move(u);
        u = std::move(u_new);
        u_new = std::move(temp_vec);

        return std::span<const double>(u);
    }

    void apply_bcs(GridPoints & buffer) const {
        if constexpr (BC_L::bc_type == NEUMANN)
            // This is a Neumann BC
            buffer[0] = buffer[1] - h * BC_L::evaluate(t);
        else
            // This is a Dirichlet BC
            buffer[0] = BC_L::evaluate(t);

        if constexpr (BC_R::bc_type == NEUMANN)
            // This is a Neumann BC
            buffer[buffer.size() - 1] = buffer[buffer.size() - 2] + h * BC_R::evaluate(t);
        else
            // This is a Dirichlet BC
            buffer[buffer.size() - 1] = BC_R::evaluate(t);
    }

private:
    long M;
    double nu, t;
    GridPoints u0;
    GridPoints u;
    GridPoints u_new;
    double h, dt;
};

#endif //HEAT_EQUATION_1D_ONE_SIDED_H
