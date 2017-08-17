
#include "parameters.h"


// NOT COMPLETED AT ALL
/* Gets the energy a configuration of positions has in a given lattice */
double get_energy(double *positions, int lattice_no, Parameters *params) {
    return 0;
}


// NOT COMPLETED AT ALL
/* Fills an array of forces with the forces on particles at given positions in a certain lattice */
void get_forces(double *forces, double *positions, int lattice_no, Parameters *params) {
    for (long d = 0; d < params->no_dimensions; d++) {
        forces[d] = 0;
    }
}


/* Implements a BAOAB step given an array of positions and forces acting on them */
void BAOAB_limit(double *current_positions, double *other_positions, double *forces, Parameters *params) {
    // Creates new array of normally distributed values
    double *new_normal_dist_arr = malloc(sizeof(double) * params->no_dimensions);
    fill_arr_normal_dist(new_normal_dist_arr, params->no_dimensions);

    for (long d = 0; d < params->no_dimensions; d++) {
        // Updates all the positions using a BAOAB step
        double change_in_position = -params->timestep * forces[d] / params->mass + \
        sqrt(0.5 * params->temp * params->timestep / params->mass) * (params->normal_dist_arr[d] + new_normal_dist_arr[d]);
        current_positions[d] += change_in_position;
        other_positions[d] += change_in_position;
    }

    // Swaps array pointers, so the array in the Parameters struct now points to the new array
    swap_arr_pointers(&params->normal_dist_arr, &new_normal_dist_arr);
    free(new_normal_dist_arr); // Frees the old array (now pointed to by new_normal_dist_arr)
}


/* Attempts a Monte-Carlo lattice switch from the current lattice to the other */
void lattice_switch(double *current_positions, double *other_positions, Parameters *params) {
    int other_lattice = (params->current_lattice + 1) % 2; // (Note current_lattice = 0 or 1

    // Calculates energies in current lattice and in the other lattice
    double current_energy = get_energy(current_positions, params->current_lattice, params);
    double other_energy   = get_energy(current_positions, other_lattice,           params);

    // Applies the energy shift to lattice 1
    if (params->current_lattice == 1) {
        current_energy += params->energy_shift;
    } else {
        other_energy += params->energy_shift;
    }

    double energy_difference = other_energy - current_energy; // Difference in energy between lattices

    // Makes the lattice switch with probability min(1, exp(-energy_difference / temp)
    if (genrand_real1() < exp(-energy_difference / params->temp)) {
        params->current_lattice = other_lattice; // Changes current_lattice
        swap_arr_pointers(&current_positions, &other_positions); // Swaps pointers of arrays for different positions
    }
}


int main(void) {
    return 0;
}