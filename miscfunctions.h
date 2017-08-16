
#ifndef BAOAB_LATTICE_SWITCH_MISCFUNCTIONS_H
#define BAOAB_LATTICE_SWITCH_MISCFUNCTIONS_H

#include <stdio.h>
#include <stdlib.h>

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



#endif //BAOAB_LATTICE_SWITCH_MISCFUNCTIONS_H
