#include <iostream>
#include <numeric>

#include <Eigen/Eigen>

enum BCType {
    DIRICHLET, NEUMANN
};

using GridPoints = Eigen::Vector<double, Eigen::Dynamic>;

// BC_L, BC_R are functors that compute boundary conditions as a function of time
template <typename BC_L, typename BC_R, typename IC>
class HeatEquation1D {
public:
    HeatEquation1D(long init_m, double init_nu, GridPoints init_u0, double init_h, double init_dt)
        : m{init_m}
        , nu{init_nu}
        , u0{std::move(init_u0)}
        , u_size{u0.size()}
        , h{init_h}
        , dt{init_dt}
    {
        if constexpr (BC_L::type == NEUMANN) {
            u_size += 1;
        }
        if constexpr (BC_R::type == NEUMANN) {
            u_size += 1;
        }
    }

private:
    long m;
    double nu;
    GridPoints u0;
    GridPoints u;
    long u_size;
    double h, dt;
};

struct BC_L {
    static constexpr BCType bc_type = DIRICHLET;
    static constexpr double evaluate(double t [[maybe_unused]]) {
        return 0.0;
    }
};

struct BC_R {
    static constexpr BCType bc_type = NEUMANN;
    static constexpr double evaluate(double t [[maybe_unused]]) {
        return 0.0;
    }
};

struct IC {
    static double evaluate(double x) {
        return cos(M_PI * x / 2.0);
    }
};

int main() {
    HeatEquation1D<BC_L, BC_R, IC> heat_equation{};

    return 0;
}
