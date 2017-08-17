
#ifndef BAOAB_LATTICE_SWITCH_MISCFUNCTIONS_H
#define BAOAB_LATTICE_SWITCH_MISCFUNCTIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mt19937ar.h"

#define PI 3.14159265358979323846264338327

/* Function to use fscanf to read doubles from a line in an input file */
void read_double(FILE *input_file, double *double_value) {
    /* Attempts to read line of input file */
    if (fscanf(input_file, "%lf", double_value) != 1) {
        printf("Failed to read parameter\n");
        exit(1);
    }
}


/* Function to use fscanf to read longs from a line in an input file */
void read_long(FILE *input_file, long *long_value) {
    /* Attempts to read line of input file */
    if (fscanf(input_file, "%ld", long_value) != 1) {
        printf("Failed to read parameter\n");
        exit(1);
    }
}


/* Function to use fscanf to read ints from a line in an input file */
void read_int(FILE *input_file, int *int_value) {
    /* Attempts to read line of input file */
    if (fscanf(input_file, "%d", int_value) != 1) {
        printf("Failed to read parameter\n");
        exit(1);
    }
}


/* Returns a random normally distributed number, mean 0, standard deviation 1 */
double box_muller_rand(){
    double r1 = genrand_real3();
    double r2 = genrand_real3();
    return sqrt(-2 * log(r1)) * cos(2 * PI * r2);
}


/* Fills an array with normally distributed values */
void fill_arr_normal_dist(double *arr, long arr_length) {
    for (long index = 0; index < arr_length; index++) {
        arr[index] = box_muller_rand();
    }
}


/* Swaps the pointers for two arrays */
void swap_arr_pointers(double **array1, double **array2) {
    double *temp;
    temp = *array1;
    *array1 = *array2;
    *array2 = temp;
}


/* Finds the minimum value of two doubles */
double min(double num1, double num2) {
    double return_val;
    if (num1 < num2) {
        return_val = num1;
    } else {
        return_val = num2;
    }

    return return_val;
}


#endif //BAOAB_LATTICE_SWITCH_MISCFUNCTIONS_H
