//
// Created by egor on 10/3/24.
//

#include "heat_equation_fwd.h"

void do_step_1st_order_bc(double ** v, double ** buffer, double c, long M) {
#pragma omp parallel for default(none) shared(M, c, v, buffer)
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
#pragma omp parallel for default(none) shared(M, c, v, buffer)
    for (long m = 1; m < M; m ++) {
        (*buffer)[m] = (*v)[m] + c * ((*v)[m + 1] - 2.0 * (*v)[m] + (*v)[m - 1]);
    }

    (*buffer)[M] = 0.0; // Dirichlet BC
    (*buffer)[0] = (*v)[0] + 2.0 * c * ((*v)[1] - (*v)[0]); // Neumann BC

    // Swap buffers
    double * temp;
    temp = *v;
    *v = *buffer;
    *buffer = temp;
}
