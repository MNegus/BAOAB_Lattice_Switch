
#include "parameters.h"
#include <time.h>

// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// THESE FUNCTIONS WILL NEED TO BE CHANGED DEPENDING ON HOW YOU ARE PASSING IN POSITIONS, ENERGIES AND FORCES

#include "example_potential.h" // Gives an example 2D potential, not needed in general

/* Gets the energy a configuration of positions has in a given lattice */
double get_unbiased_energy(double *positions, int current_lattice) {
    double return_val;
    if (current_lattice == 0) {
        return_val = potential_1(positions);
    } else {
        return_val = potential_2(positions);
    }
    return return_val;
}


/* Fills an array of forces with the forces on particles at given positions in a certain lattice */
void get_unbiased_forces(double *forces, double *positions, int current_lattice) {
    if (current_lattice == 0) {
        fill_potential_1_deriv(forces, positions);
    } else {
        fill_potential_2_deriv(forces, positions);
    }
}


/* Fills an array of positions with the required initial positions of each lattice */
void get_initial_positions(double *lattice1_positions, double *lattice2_positions) {
    lattice1_positions[0] = 0;
    lattice1_positions[1] = 0;

    lattice2_positions[0] = 100;
    lattice2_positions[1] = 0;
}
// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/* Difference in the energy between lattice 1 and lattice 2 */
double energy_diff_one_to_two(double *lattice1_positions, double *lattice2_positions, Parameters *params) {
    // Returns the difference in the unbiased energy in lattice 1 and lattice 2, minus the shift from lattice 2
    return get_unbiased_energy(lattice1_positions, 0) - get_unbiased_energy(lattice2_positions, 1) -
           params->energy_shift;
}


/* Gets the value of the biasing function for a certain configuration */
double biasing_function(double *lattice1_positions, double *lattice2_positions, Parameters *params) {
    double bias_func_val = 0; // Value of the biasing function at the current position

    // Energy difference between lattice 1 and lattice 2
    double energy_diff = energy_diff_one_to_two(lattice1_positions, lattice2_positions, params);

    // Adds contributions from all the Gaussians
    for (long gauss_no = 0; gauss_no < params->no_gaussians; gauss_no++) {
        bias_func_val += gaussian(energy_diff, params->gauss_width, params->gauss_positions[gauss_no],
                                  params->gauss_heights[gauss_no]);
    }

    return bias_func_val;
}


/* Get the energy of the current lattice after biasing has been applied by the Gaussians */
double get_biased_energy(int current_lattice, double *lattice1_positions, double *lattice2_positions,
                         Parameters *params) {
    double return_val = params->temp * biasing_function(lattice1_positions, lattice2_positions, params);

    if (current_lattice == 0) {
        return_val += get_unbiased_energy(lattice1_positions, current_lattice);
    } else {
        return_val += get_unbiased_energy(lattice2_positions, current_lattice) + params->energy_shift;
    }

    return return_val;
}


void get_biased_forces(int current_lattice, double *bias_forces, double *lattice1_positions, double *lattice2_positions,
                       Parameters *params) {

    fill_double_arr_zeros(bias_forces, params->no_dimensions); // Ensures the biased forces array is empty

    // Finds forces on each lattice
    double *lattice1_forces = malloc(sizeof(double) * params->no_dimensions);
    get_unbiased_forces(lattice1_forces, lattice1_positions, 0);

    double *lattice2_forces = malloc(sizeof(double) * params->no_dimensions);
    get_unbiased_forces(lattice2_forces, lattice2_positions, 1);

    // Determines which array points to forces in the current lattice
    double *current_forces;
    if (current_lattice == 0) {
        current_forces = lattice1_forces;
    } else {
        current_forces = lattice2_forces;
    }

    // Energy difference between lattice 1 and lattice 2
    double energy_diff = energy_diff_one_to_two(lattice1_positions, lattice2_positions, params);

    // Loops over all dimensions and calculates the biased force
    for (long d = 0; d < params->no_dimensions; d++) {
        // Adds the contributions from the derivative of the Gaussians
        for (long gauss_no = 0; gauss_no < params->no_gaussians; gauss_no++) {
            bias_forces[d] += gaussian_deriv(energy_diff, params->gauss_width,
                                             params->gauss_positions[gauss_no], params->gauss_heights[gauss_no]);
        }
        // Biased force is difference in forces multiplied by sum of Gaussians multiplied by temp, plus the unbiased force
        bias_forces[d] = params->temp * (lattice1_forces[d] - lattice2_forces[d]) * bias_forces[d] + current_forces[d];
    }

    free(lattice1_forces);
    free(lattice2_forces);
}


/* Implements a BAOAB step given an array of positions and forces acting on them */
void BAOAB_limit(double *lattice1_positions, double *lattice2_positions, double *forces, Parameters *params) {
    // Creates new array of normally distributed values
    double *new_normal_dist_arr = malloc(sizeof(double) * params->no_dimensions);
    fill_arr_normal_dist(new_normal_dist_arr, params->no_dimensions);

    for (long d = 0; d < params->no_dimensions; d++) {
        // Updates all the positions using a BAOAB step
        double change_in_position = -params->timestep * forces[d] / params->mass + \
        sqrt(0.5 * params->temp * params->timestep / params->mass) * (params->normal_dist_arr[d] +
                                                                      new_normal_dist_arr[d]);
        lattice1_positions[d] += change_in_position;
        lattice2_positions[d] += change_in_position;
    }

    // Swaps array pointers, so the array in the Parameters struct now points to the new array
    swap_arr_pointers(&params->normal_dist_arr, &new_normal_dist_arr);
    free(new_normal_dist_arr); // Frees the old array (now pointed to by new_normal_dist_arr)
}


/* Attempts a Monte-Carlo lattice switch from the current lattice to the other */
void lattice_switch(double *lattice1_positions, double *lattice2_positions, Parameters *params) {
    int other_lattice = (params->current_lattice + 1) % 2;

    double current_biased_energy = get_biased_energy(params->current_lattice, lattice1_positions, lattice2_positions,
                                                     params);
    double other_biased_energy = get_biased_energy(other_lattice, lattice1_positions, lattice2_positions, params);

    // Makes the lattice switch with probability min(1, exp(-relative_energy_diff / temp)
    if (genrand_real1() < exp(-(other_biased_energy - current_biased_energy) / params->temp)) {
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

        get_unbiased_forces(forces, lattice1_positions, 0); // Calculates forces from the first lattice
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
        get_unbiased_forces(forces, lattice2_positions, 1);
        BAOAB_limit(lattice2_positions, lattice1_positions, forces, params);
    }

    swap_arr_pointers(&initial_lattice1_positions, &lattice1_positions);
    free(initial_lattice1_positions);
    swap_arr_pointers(&initial_lattice2_positions, &lattice2_positions);
    free(initial_lattice2_positions);

    free(forces);
}


/* Implements the biasing, so that we have arrays of heights and positions of Gaussians in energy-difference space, placed
 * such that all points in energy-difference space are equally likely to be visited */
void implement_biasing(double *lattice1_positions, double *lattice2_positions, Parameters *params) {

    // Variable to store the unbiased energy difference between lattice 1 and lattice 2
    double energy_diff = energy_diff_one_to_two(lattice1_positions, lattice2_positions, params);

    // Creates an empty histogram to record the energy differences we have visited
    long *energy_diff_histogram = malloc(sizeof(double) * params->no_bins);
    fill_long_arr_zeros(energy_diff_histogram, params->no_bins);

    // Defines variable f, used to alter the height of the Gaussians
    double f = exp(params->initial_gauss_height / params->temp); // Initial value for f
    double min_f = exp(params->min_gauss_height / params->temp); // Smallest value f is allowed to take

    double *bias_forces = malloc(
            sizeof(double) * params->no_dimensions); // Array for the bias forces on a particular lattice


    // Decreases f until its at its minimum. At that point the Gaussians are sufficient for our biasing
    while (f > min_f) {
        int flat_histogram = 0; // Indicates if our histogram is flat (0 being not flat, 1 being flat)

        // Adds Gaussians until our histogram becomes sufficiently flat
        while (flat_histogram == 0) {



            // Performs a set number of BAOAB steps
            for (long stepno = 0; stepno < params->steps_between_gaussians; stepno++) {

                // Adds to the histogram with the current energy difference
                add_to_histogram(energy_diff, energy_diff_histogram, params->no_bins, params->minimum_energy_diff,
                                 params->maximum_energy_diff, params->bin_width);

                // Attempts a lattice switch with set probability
                if (genrand_real1() < params->prob_switch_attempt) {
                    lattice_switch(lattice1_positions, lattice2_positions, params);
                } else {
                    // Else, gets the biased forces on the current lattice and performs a BAOAB step
                    get_biased_forces(params->current_lattice, bias_forces, lattice1_positions, lattice2_positions,
                                      params);
                    BAOAB_limit(lattice1_positions, lattice2_positions, bias_forces, params);

                    // Updates the energy difference (only done on BAOAB as the lattice switch does not change this)
                    energy_diff = energy_diff_one_to_two(lattice1_positions, lattice2_positions, params);
                }
            }

            // Places a Gaussian where we currently are
            params->gauss_positions[params->no_gaussians] = energy_diff;
            params->gauss_heights[params->no_gaussians] = params->temp * log(f);
            params->no_gaussians++;

            // If the number of Gaussians now exceeds the length of the height and position arrays, extend their length
            if (params->no_gaussians >= params->gauss_arrs_length) {
                long new_arrs_length = params->gauss_arrs_length + params->init_gauss_arrs_length;
                params->gauss_positions = extend_arr(params->gauss_positions, params->gauss_arrs_length,
                                                     new_arrs_length);
                params->gauss_heights = extend_arr(params->gauss_heights, params->gauss_arrs_length, new_arrs_length);
                params->gauss_arrs_length += params->init_gauss_arrs_length;
            }

            // Determines if the histogram is now flat or not
            flat_histogram = is_histogram_flat(energy_diff_histogram, params->no_bins, params->flatness_threshold);
        }

        f = sqrt(f); // Decrease f
        fill_long_arr_zeros(energy_diff_histogram, params->no_bins); // Reset the histogram
    }

    free(bias_forces);
    free(energy_diff_histogram);
}


int main(int argc, char **argv) {
    char *input_filename = argv[1]; // Name of input file for parameters (see parameters.h for details)
    char *output_filename = argv[2]; // Output file for free energy differences at every timestep

    // /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Bad way of generating the seed when running lots of simulations at once
    unsigned long seed = (unsigned long) time(NULL);
    // /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    init_genrand(seed);

    Parameters params;
    store_parameters(&params, input_filename);

    double *initial_lattice1_positions = malloc(sizeof(double) * params.no_dimensions);
    double *initial_lattice2_positions = malloc(sizeof(double) * params.no_dimensions);

    get_initial_positions(initial_lattice1_positions, initial_lattice2_positions);

    params.energy_shift = get_unbiased_energy(initial_lattice1_positions, 0) - \
        get_unbiased_energy(initial_lattice2_positions, 1);

    double *lattice1_positions = malloc(sizeof(double) * params.no_dimensions);
    double *lattice2_positions = malloc(sizeof(double) * params.no_dimensions);

    copy_arr(initial_lattice1_positions, lattice1_positions, params.no_dimensions);
    copy_arr(initial_lattice2_positions, lattice2_positions, params.no_dimensions);

    implement_biasing(lattice1_positions, lattice2_positions, &params);

    copy_arr(initial_lattice1_positions, lattice1_positions, params.no_dimensions);
    copy_arr(initial_lattice2_positions, lattice2_positions, params.no_dimensions);
    params.current_lattice = params.start_lattice;

    double biased_lattice1_count = 0;
    double biased_lattice2_count = 0;

    double free_energy_difference;

    double *bias_forces = malloc(sizeof(double) * params.no_dimensions);

    FILE *output_file = fopen(output_filename, "w");

    for (long step_no = 0; step_no < params.tot_timesteps; step_no++) {
        double bias_func_val = biasing_function(lattice1_positions, lattice2_positions, &params);
        if (params.current_lattice == 0) {
            biased_lattice1_count += exp(bias_func_val);
        } else {
            biased_lattice2_count += exp(bias_func_val);
        }

        free_energy_difference =
                -params.temp * log(biased_lattice1_count / biased_lattice2_count) + params.energy_shift;
        fprintf(output_file, "%ld, %lf\n", step_no, free_energy_difference);

        if (genrand_real1() < params.prob_switch_attempt) {
            lattice_switch(lattice1_positions, lattice2_positions, &params);
        } else {
            get_biased_forces(params.current_lattice, bias_forces, lattice1_positions, lattice2_positions, &params);
            BAOAB_limit(lattice1_positions, lattice2_positions, bias_forces, &params);
        }

    }

    fclose(output_file);

    free(initial_lattice1_positions);
    free(initial_lattice2_positions);
    free(lattice1_positions);
    free(lattice2_positions);
    free(bias_forces);
    free_parameters(&params);

    return 0;
}