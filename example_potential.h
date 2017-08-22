//
// Created by michael on 17/08/17.
//


#ifndef BAOAB_LATTICE_SWITCH_EXAMPLE_POTENTIAL_H
#define BAOAB_LATTICE_SWITCH_EXAMPLE_POTENTIAL_H

#include <float.h>


/* Example potential for lattice 1*/
double potential_1(double *positions) {
    return 10 * positions[0] * positions[0] + positions[1] * positions[1];
}


/* Example potential for lattice 2 */
double potential_2(double *positions) {
    return (positions[0] - 100) * (positions[0] - 100) + 10 * positions[1] * positions[1] - 24;
}


/* Fills an array with the derivatives for the lattice 1 potential */
void fill_potential_1_deriv(double *deriv_arr, double *positions) {
    deriv_arr[0] = 20 * positions[0];
    deriv_arr[1] = 2 * positions[1];
}


/* Fills an array with the derivatives for the lattice 2 potential */
void fill_potential_2_deriv(double *deriv_arr, double *positions) {
    deriv_arr[0] = 2 * (positions[0] - 100);
    deriv_arr[1] = 20 * positions[1];
}

#endif //BAOAB_LATTICE_SWITCH_EXAMPLE_POTENTIAL_H
