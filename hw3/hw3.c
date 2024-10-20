//
// Created by egor on 10/19/24.
//

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

double p_i_next(double p_i, double a_i, double b_i, double c_i) {
    return -c_i / (a_i * p_i + b_i);
}

double q_i_next(double q_i, double p_i, double a_i, double b_i, double d_i) {
    return (d_i - a_i * q_i) / (q_i * p_i + b_i);
}


int main() {

    const double a = 1.0;
    const double b = 0.0;
    const double c = 1.0;
    const double d = 0.0;

    const long M = 10;

    const double beta_0 = 0.0;

    return 0;
}
