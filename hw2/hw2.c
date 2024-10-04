//
// Created by egor on 10/2/24.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "heat_equation_fwd.h"
#include "../csv_writer.h"

void part1() {
    const long M = 10;
    const double h = 1.0 / (double) M;
    const double dt = 0.004;
    const double nu = 1.0;
    const double c = nu * dt / h / h;
    const double time_points[] = {0.06, 0.1, 0.9};

    // Allocate memory for the solution vector and the temporary buffer
    double * v, * buffer;
    v = (double *) malloc(sizeof(double) * (M + 1));
    buffer = (double *) malloc(sizeof(double) * (M + 1));

    // Populate the initial condition
    for (long m = 0; m < M + 1; m ++) {
        double xm = h * (double) m;
        v[m] = cos(M_PI * xm / 2.0);
    }

    for (long n = 1; n < 1000000; n ++) {
        do_step_2nd_order_bc(&v, &buffer, c, M);
    }

    // Output the results
    for (long m = 0; m < M + 1; m ++) {
        double xm = h * (double) m;
        printf("%.4f\t%.4f\n", xm, v[m]);
    }

    // Free the buffers
    free(v);
    free(buffer);
}

double compute_l2_norm_of_error(double * v, double t, double nu, long M) {
    double l2_norm = 0.0;
    for (long m = 0; m < M + 1; m ++) {
        double xm = (double) m / (double) M;
        double u = cos(M_PI * xm / 2.0) * exp(-nu * M_PI * M_PI / 4.0 * t);
        l2_norm += (v[m] - u) * (v[m] - u);
    }
    return sqrt(l2_norm);
}

void part2() {
    const long M_values[] = {8, 16, 32, 64, 128, 256, 512/*, 1024, 2048*/};
    const double dt_over_h2 = 0.004 / 0.01;
    const double nu = 1.0;
    const double t_tot = 2.0;

    const char * col_names[] = {"h", "log(h)", "o1", "log(o1)", "o2", "log(o2)"};
    csv_writer_t csv_writer = open_csv_writer("hw2_part2.csv", col_names, 6);

    for (long i = 0; i < sizeof(M_values) / sizeof(M_values[0]); i ++) {
        long M = M_values[i];
        double h = 1.0 / (double) M;
        double dt = dt_over_h2 * h * h;
        double c = nu * dt / h / h;
        long N = (long) (t_tot / dt);

        double l2_norm_1o, l2_norm_2o;

        // Allocate memory for the solution vector and the temporary buffer
        double * v, * buffer;
        v = (double *) malloc(sizeof(double) * (M + 1));
        buffer = (double *) malloc(sizeof(double) * (M + 1));

        // Populate the initial condition
        for (long m = 0; m < M + 1; m ++) {
            double xm = h * (double) m;
            v[m] = cos(M_PI * xm / 2.0);
        }

        for (long n = 1; n < N; n ++) {
            do_step_1st_order_bc(&v, &buffer, c, M);
        }

        l2_norm_1o = compute_l2_norm_of_error(v, t_tot, nu, M);

        // Populate the initial condition
        for (long m = 0; m < M + 1; m ++) {
            double xm = h * (double) m;
            v[m] = cos(M_PI * xm / 2.0);
        }

        for (long n = 1; n < N; n ++) {
            do_step_2nd_order_bc(&v, &buffer, c, M);
        }

        l2_norm_2o = compute_l2_norm_of_error(v, t_tot, nu, M);

        printf("%lf\t%le\t%le\n", h, l2_norm_1o, l2_norm_2o);

        double csv_line[] = {h, log2(h), l2_norm_1o, log2(l2_norm_1o), l2_norm_2o, log2(l2_norm_2o)};
        append_csv_line(csv_writer, csv_line);

        // Free the buffers
        free(v);
        free(buffer);
    }

    close_csv_writer(csv_writer);
}

int main() {
    part2();

    return 0;
}
