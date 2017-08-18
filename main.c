
#include "parameters.h"
#include "example_potential.h"
#include <float.h>
#include <time.h>


// NOT COMPLETED AT ALL
/* Gets the energy a configuration of positions has in a given lattice */
double get_unbiased_energy(double *positions, Parameters *params) {
    return potential(positions);
}


/* Difference in the energy between lattice 1 and lattice 2 */
double energy_diff_one_to_two(double *lattice1_positions, double *lattice2_positions, Parameters *params) {
    return get_unbiased_energy(lattice1_positions, params) - get_unbiased_energy(lattice2_positions, params) - params->energy_shift;
}


/* Relative energy difference between the current lattice and the other one */
double relative_energy_diff(double *current_positions, double *other_positions, Parameters *params) {
    double energy_diff;
    if (params->current_lattice == 0) {
        energy_diff = energy_diff_one_to_two(current_positions, other_positions, params);
    } else {
        energy_diff = energy_diff_one_to_two(other_positions, current_positions, params);
    }

    return energy_diff;
}


/* Get the energy of the current lattice after biasing has been applied by the Gaussians */
double biased_energy(double *lattice1_positions, double *lattice2_positions,
                     double *gauss_heights, double *gauss_positions, long no_gaussians, Parameters *params) {
    double *current_positions, *other_positions;
    if (params->current_lattice == 0) {
        current_positions = lattice1_positions;
        other_positions   = lattice2_positions;
    } else {
        current_positions = lattice2_positions;
        other_positions   = lattice1_positions;
    }

    double bias_func_val = 0;

    double energy_diff = relative_energy_diff(current_positions, other_positions, params);

    for (long gauss_no = 0; gauss_no < no_gaussians; gauss_no++) {
        bias_func_val += gaussian(energy_diff, params->gauss_width, gauss_positions[gauss_no], gauss_heights[gauss_no]);
    }

    return params->temp * bias_func_val + get_unbiased_energy(current_positions, params);
}


// NOT COMPLETED AT ALL
/* Fills an array of forces with the forces on particles at given positions in a certain lattice */
void get_unbiased_forces(double *forces, double *positions, Parameters *params) {
    potential_deriv(forces, positions);
}


void biased_forces(double *forces, double *lattice1_positions, double *lattice2_positions,
                   double *gauss_heights, double *gauss_positions, long no_gaussians, Parameters *params) {
    
}


/* Implements a BAOAB step given an array of positions and forces acting on them */
void BAOAB_limit(double *lattice1_positions, double *lattice2_positions, double *forces, Parameters *params) {
    // Creates new array of normally distributed values
    double *new_normal_dist_arr = malloc(sizeof(double) * params->no_dimensions);
    fill_arr_normal_dist(new_normal_dist_arr, params->no_dimensions);

    for (long d = 0; d < params->no_dimensions; d++) {
        // Updates all the positions using a BAOAB step
        double change_in_position = -params->timestep * forces[d] / params->mass + \
        sqrt(0.5 * params->temp * params->timestep / params->mass) * (params->normal_dist_arr[d] + new_normal_dist_arr[d]);
        lattice1_positions[d] += change_in_position;
        lattice2_positions[d] += change_in_position;
    }

    // Swaps array pointers, so the array in the Parameters struct now points to the new array
    swap_arr_pointers(&params->normal_dist_arr, &new_normal_dist_arr);
    free(new_normal_dist_arr); // Frees the old array (now pointed to by new_normal_dist_arr)
}


/* Attempts a Monte-Carlo lattice switch from the current lattice to the other */
void lattice_switch(double *lattice1_positions, double *lattice2_positions, Parameters *params) {
    double energy_diff;
    if (params->current_lattice == 0) {
        energy_diff = energy_diff_one_to_two(lattice1_positions, lattice2_positions, params);
    } else {
        energy_diff = -energy_diff_one_to_two(lattice1_positions, lattice2_positions, params);
    }

    // Makes the lattice switch with probability min(1, exp(-energy_diff_one_to_two / temp)
    if (genrand_real1() < exp(-energy_diff / params->temp)) {
        params->current_lattice = (params->current_lattice + 1) % 2;  // Changes current_lattice
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
        double current_energy_difference = energy_diff_one_to_two(lattice1_positions, lattice2_positions, params);

        // Checks if current energy difference is new minimum or maximum
        if (current_energy_difference < params->minimum_energy_diff) {
            params->minimum_energy_diff = current_energy_difference;
        }
        if (current_energy_difference > params->maximum_energy_diff) {
            params->maximum_energy_diff = current_energy_difference;
        }

        get_unbiased_forces(forces, lattice1_positions, params); // Calculates forces from the first lattice
        BAOAB_limit(lattice1_positions, lattice2_positions, forces, params);
    }

    // Initialised arrays back to their original positions
    copy_arr(initial_lattice1_positions, lattice1_positions, params->no_dimensions);
    copy_arr(initial_lattice2_positions, lattice2_positions, params->no_dimensions);


    // Runs BAOAB stepping simulation on the second lattice
    for (long step_no = 0; step_no < params->tot_timesteps; step_no++) {
        double current_energy_difference = -energy_diff_one_to_two(lattice1_positions, lattice2_positions, params);
        if (current_energy_difference < params->minimum_energy_diff) {
            params->minimum_energy_diff = current_energy_difference;
        }
        if (current_energy_difference > params->maximum_energy_diff) {
            params->maximum_energy_diff = current_energy_difference;
        }
        get_unbiased_forces(forces, lattice2_positions, params);
        BAOAB_limit(lattice2_positions, lattice1_positions, forces, params);
    }

    swap_arr_pointers(&initial_lattice1_positions, &lattice1_positions);
    free(initial_lattice1_positions);
    swap_arr_pointers(&initial_lattice2_positions, &lattice2_positions);
    free(initial_lattice2_positions);

    free(forces);
}


void implement_biasing(double *gauss_heights, double *gauss_positions,
                       double *lattice1_positions, double *lattice2_positions, Parameters *params) {
    // Creates an empty histogram to record the energy differences we have visited
    long *energy_diff_histogram = malloc(sizeof(double) * params->no_bins);
    fill_arr_zeros(energy_diff_histogram, params->no_bins);

    // Defines variable f, used to alter the height of the Gaussians
    double f = exp(params->initial_gauss_height / params->temp); // Initial value for f
    double min_f = exp(params->min_gauss_height / params->temp); // Smallest value f is allowed to take

    long no_gaussians = 0; // Number of Gaussians that have currently been placed

    while (f > min_f) {
        int flat_histogram = 0; // Indicates if our histogram is flat (0 being not flat, 1 being flat)

        while (flat_histogram == 0) {

        }
    }




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

    free(initial_energies);
    free(lattice1_positions);
    free(lattice2_positions);

    return 0;
}