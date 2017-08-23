
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
double box_muller_rand() {
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


/* Copies the contents of the array original_arr into the array arr_to_be_filled */
void copy_arr(double *original_arr, double *arr_to_be_filled, long arr_length) {
    for (long index = 0; index < arr_length; index++) {
        arr_to_be_filled[index] = original_arr[index];
    }

}


/* "Extends" the length of an array by creating a new, longer array, and copying the contents from the old one into
 * it, and then returning a pointer to the new one */
double *extend_arr(double *current_arr, long current_arr_length, long new_arr_length) {
    double *new_arr = malloc(sizeof(double) * new_arr_length); // New array
    copy_arr(current_arr, new_arr, current_arr_length); // Copies contents from old to new
    free(current_arr); // Frees the old array
    return new_arr; // Return pointer to new array
}

/* Fills an array of longs with zeros */
void fill_long_arr_zeros(long *arr, long arr_length) {
    for (long index = 0; index < arr_length; index++) {
        arr[index] = 0;
    }
}


/* Fills an array of doubles with zeros */
void fill_double_arr_zeros(double *arr, long arr_length) {
    for (long index = 0; index < arr_length; index++) {
        arr[index] = 0;
    }
}


/* Gaussian function with given width, height and mean */
double gaussian(double position, double width, double mean, double height) {
    return height * exp(-(position - mean) * (position - mean) / width);
}


/* Derivative of a gaussian with given width, height and mean */
double gaussian_deriv(double position, double width, double mean, double height) {
    return -2 * (position - mean) * gaussian(position, width, mean, height) / width;
}


/* Gets the number of the bin which a given position belongs in for a histogram */
long get_bin_no(double position, long no_bins, double min_pos, double max_pos, double bin_width) {
    long bin_no; // Variable to be returned
    if (position < min_pos) {
        // If the position is to the left of the first bin, then it belongs in bin 0
        bin_no = -1;
    } else if (position > max_pos) {
        // If the position is to the right of the last bin, then it belongs in the last bin
        bin_no = -1;
    } else {
        // Loops over all bins until we find which it lies in, and then breaks the loop
        for (long j = 1; j <= no_bins; j++) {
            if (position < min_pos + j * bin_width) {
                bin_no = j - 1;
                break;
            }
        }
    }

    return bin_no;
}


/* Given a position, adds to the bin it belongs in in a histogram */
void add_to_histogram(double position, long *hist, long no_bins, double min_pos, double max_pos, double bin_width) {
    long bin_no = get_bin_no(position, no_bins, min_pos, max_pos, bin_width);
    if (bin_no != -1) hist[bin_no]++;
}


/* Determines if a histogram is sufficiently "flat", meaning that all bins have a value greater than or
 * equal to flatness_threshold * mean */
int is_histogram_flat(long *hist, long no_bins, double flatness_threshold) {
    // Calculates mean of the bin values
    double mean = 0;
    for (long bin_no = 0; bin_no < no_bins; bin_no++) {
        mean += hist[bin_no];
    }
    mean /= no_bins;

    int histogram_flat = 1; // If 1, then histogram flat. If 0, then it is not

    /* Loops through all bins. If any bin has a value less than the required for flatness, we set histogram_flat
     * to be zero and break out of the loop */
    for (long bin_no = 0; bin_no < no_bins; bin_no++) {
        if (hist[bin_no] < flatness_threshold * mean) {
            histogram_flat = 0;
            break;
        }
    }

    return histogram_flat;
}

#endif //BAOAB_LATTICE_SWITCH_MISCFUNCTIONS_H
