//
// Created by egor on 10/19/24.
//

#include <assert.h>
#include <stdlib.h>

void do_step(double aa, double bb, double cc, double a, double lambda, double * p, double * q, double * v, long M) {

    assert(p);
    assert(q);
    assert(v);

    p[0] = 0.0;
    q[1] = 0.0;  // Data(time)

    for (long m = 1; m < M - 1; m ++) {
        double dd = v[m] - a*lambda*(v[m+1] - v[m-1]) / 2.0;
        double denom = aa * p[m] + bb;
        p[m+1] = -cc / denom;
        q[m+1] = (dd - q[m] * aa) / denom;
    }

    // Apply boundary conditions at the last point
    v[M - 1] = v[M - 2];

    // Compute all interior values
    for (long m = M - 1; m > 0; m --) {
        v[m] = p[m+1] * v[m+1] + q[m+1];
    }
}

int main() {
    const long M = 10;

    const double a = 1.0;
    const double lambda = 1.0;
    const double aa = -a*lambda/4.0;
    const double bb = 1.0;
    const double cc = -aa;

    // Set the first elements of the p and q arrays
    double * p = (double *) malloc(sizeof(double) * (M - 1));
    double * q = (double *) malloc(sizeof(double) * (M - 1));
    double * v = (double *) malloc(sizeof(double) * M);

    // Initialize v to the IC


    do_step(aa, bb, cc, a, lambda, p, q, v, M);

    free(p);
    free(q);
    free(v);

    return 0;
}
