//
// Created by michael on 17/08/17.
//

#include "parameters.h"

#ifndef BAOAB_LATTICE_SWITCH_EXAMPLE_POTENTIAL_H
#define BAOAB_LATTICE_SWITCH_EXAMPLE_POTENTIAL_H


/* Example two dimensional potential */
double potential(double *positions) {

    double return_val;

    if (positions[0] < 10) {
        return_val = 10 * positions[0] * positions[0] + positions[1] * positions[1];
    } else if (positions[0] < 68) {
        return_val = 0;
    } else {
        return_val = (positions[0] - 100) * (positions[0] - 100) + 10 * positions[1] * positions[1] - 24;
    }

    return return_val;
}

/* Derivative of the example potential */
void potential_deriv(double *deriv_arr, double *positions) {


    if (positions[0] < 10) {
        deriv_arr[0] = 20 * positions[0];
        deriv_arr[1] = 2 * positions[1];
    } else if (positions[0] < 68) {
        deriv_arr[0] = 0;
        deriv_arr[1] = 0;
    } else {
        deriv_arr[0] = 2 * (positions[0] - 100);
        deriv_arr[1] = 20 * positions[1];
    }
}


#endif //BAOAB_LATTICE_SWITCH_EXAMPLE_POTENTIAL_H
