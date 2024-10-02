//
// Created by egor on 10/2/24.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void do_step_1st_order_bc(double ** v, double ** buffer, double c, long M) {
    #pragma omp parallel for
    for (long m = 1; m < M; m ++) {
        (*buffer)[m] = (*v)[m] + c * ((*v)[m + 1] - 2.0 * (*v)[m] + (*v)[m - 1]);
    }
    (*buffer)[M] = 0.0; // Dirichlet BC
    (*buffer)[0] = (*buffer)[1]; // Neumann BC

    // Swap buffers
    double * temp;
    temp = *v;
    *v = *buffer;
    *buffer = temp;
}

void do_step_2nd_order_bc(double ** v, double ** buffer, double c, long M) {
#pragma omp parallel for
    for (long m = 1; m < M; m ++) {
        (*buffer)[m] = (*v)[m] + c * ((*v)[m + 1] - 2.0 * (*v)[m] + (*v)[m - 1]);
    }
    (*buffer)[M] = 0.0; // Dirichlet BC
    (*buffer)[0] = (*buffer)[1]; // Neumann BC

    // Swap buffers
    double * temp;
    temp = *v;
    *v = *buffer;
    *buffer = temp;
}

int main() {
    const long M = 10000;
    const double h = 1.0 / (double) M;
    const double dt = 0.000000004;
    const double nu = 1.0;
    const double c = nu * dt / h / h;

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
        do_step_1st_order_bc(&v, &buffer, c, M);
    }

    // Output the results
    for (long m = 0; m < M + 1; m ++) {
        double xm = h * (double) m;
        printf("%.4f\t%.4f\n", xm, v[m]);
    }

    // Free the buffers
    free(v);
    free(buffer);

    return 0;
}
