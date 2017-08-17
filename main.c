
#include "parameters.h"
#include "example_potential.h"
#include <float.h>
#include <time.h>


// NOT COMPLETED AT ALL
/* Gets the energy a configuration of positions has in a given lattice */
double get_energy(double *positions, Parameters *params) {
    return potential(positions);
}


/* Difference in the energy between the positions in the current lattice and in the other lattice */
double energy_difference(double *current_positions, double *other_positions, Parameters *params) {
    double return_val = get_energy(other_positions, params) - get_energy(current_positions, params);
    if (params->current_lattice == 0) {
        return_val -= params->energy_shift;
    } else {
        return_val += params->energy_shift;
    }

    return return_val;
}

// NOT COMPLETED AT ALL
/* Fills an array of forces with the forces on particles at given positions in a certain lattice */
void get_forces(double *forces, double *positions, Parameters *params) {
    potential_deriv(forces, positions);
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
    // Makes the lattice switch with probability min(1, exp(-energy_difference / temp)
    if (genrand_real1() < exp(-energy_difference(current_positions, other_positions, params) / params->temp)) {
        params->current_lattice = (params->current_lattice + 1) % 2;  // Changes current_lattice
        swap_arr_pointers(&current_positions, &other_positions); // Swaps pointers of arrays for different positions
    }
}


/* Calculates the minimum and maximum energy differences that will be encountered during the simulations */
void calculate_energy_diff_range(double *lattice1_positions, double *lattice2_positions, Parameters *params) {
    // Saves the initial configurations into arrays
    double *initial_lattice1_positions = malloc(sizeof(double) * params->no_dimensions);
    double *initial_lattice2_positions = malloc(sizeof(double) * params->no_dimensions);
    copy_arr(lattice1_positions, initial_lattice1_positions, params->no_dimensions);
    copy_arr(lattice2_positions, initial_lattice2_positions, params->no_dimensions);

    // Initialised maximum and minimum energy differences
    params->minimum_energy_diff = DBL_MAX;
    params->maximum_energy_diff = -DBL_MAX;

    double *forces = malloc(sizeof(double) * params->no_dimensions); // Array for storing values of forces

    // Runs BAOAB stepping simulation on the first lattice
    for (long step_no = 0; step_no < params->tot_timesteps; step_no++) {
        // Energy difference between lattice 1 and lattice 2
        double current_energy_difference = energy_difference(lattice2_positions, lattice1_positions, params);

        // Checks if current energy difference is new minimum or maximum
        if (current_energy_difference < params->minimum_energy_diff) {
            params->minimum_energy_diff = current_energy_difference;
        }
        if (current_energy_difference > params->maximum_energy_diff) {
            params->maximum_energy_diff = current_energy_difference;
        }

        get_forces(forces, lattice1_positions, params); // Calculates forces from the first lattice
        BAOAB_limit(lattice1_positions, lattice2_positions, forces, params);
    }

    // Initialised arrays back to their original positions
    copy_arr(initial_lattice1_positions, lattice1_positions, params->no_dimensions);
    copy_arr(initial_lattice2_positions, lattice2_positions, params->no_dimensions);


    // Runs BAOAB stepping simulation on the second lattice
    for (long step_no = 0; step_no < params->tot_timesteps; step_no++) {
        double current_energy_difference = energy_difference(lattice2_positions, lattice1_positions, params);
        if (current_energy_difference < params->minimum_energy_diff) {
            params->minimum_energy_diff = current_energy_difference;
        }
        if (current_energy_difference > params->maximum_energy_diff) {
            params->maximum_energy_diff = current_energy_difference;
        }
        get_forces(forces, lattice2_positions, params);
        BAOAB_limit(lattice2_positions, lattice1_positions, forces, params);
    }

    swap_arr_pointers(&initial_lattice1_positions, &lattice1_positions);
    free(initial_lattice1_positions);
    swap_arr_pointers(&initial_lattice2_positions, &lattice2_positions);
    free(initial_lattice2_positions);

    free(forces);
}


int main(int argc, char **argv) {
    char *input_filename = argv[1];

    long seed = (unsigned long) time(NULL);
    init_genrand(seed);

    double *initial_energies = malloc(sizeof(double) * 2);
    double *lattice1_positions = malloc(sizeof(double) * 2);
    double *lattice2_positions = malloc(sizeof(double) * 2);

    lattice1_positions[0] = 0;
    lattice1_positions[1] = 0;

    lattice2_positions[0] = 100;
    lattice2_positions[1] = 0;

    initial_energies[0] = potential(lattice1_positions);
    initial_energies[1] = potential(lattice2_positions);

    Parameters params;

    store_parameters(&params, input_filename, initial_energies);

    calculate_energy_diff_range(lattice1_positions, lattice2_positions, &params);

    printf("%lf %lf\n", params.minimum_energy_diff, params.maximum_energy_diff);



    return 0;
}