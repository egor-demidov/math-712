//
// Created by egor on 10/3/24.
//

#ifndef HEAT_EQUATION_FWD_H
#define HEAT_EQUATION_FWD_H

void do_step_1st_order_bc(double ** v, double ** buffer, double c, long M);

void do_step_2nd_order_bc(double ** v, double ** buffer, double c, long M);

#endif //HEAT_EQUATION_FWD_H
