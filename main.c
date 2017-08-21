
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
    // Returns the difference in the unbiased energy in lattice 1 and lattice 2, minus the shift from lattice 2
    return get_unbiased_energy(lattice1_positions, params) - get_unbiased_energy(lattice2_positions, params) - params->energy_shift;
}


/* Relative energy difference between the current lattice and the other one */
double relative_energy_diff(int current_lattice, double *current_positions, double *other_positions, Parameters *params) {
    double energy_diff;

    // Determines which positions array belongs in which lattice, and determines the difference by calling energy_diff_one_to_two
    if (current_lattice == 0) {
        energy_diff = -energy_diff_one_to_two(current_positions, other_positions, params);
    } else {
        energy_diff = energy_diff_one_to_two(other_positions, current_positions, params);
    }

    return energy_diff;
}


/* Get the energy of the current lattice after biasing has been applied by the Gaussians */
double get_biased_energy(int current_lattice, double *lattice1_positions, double *lattice2_positions,
                         Parameters *params) {
    // Determines which positions array points to the lattice we are currently in
    double *current_positions;
    if (current_lattice == 0) {
        current_positions = lattice1_positions;
    } else {
        current_positions = lattice2_positions;
    }

    double bias_func_val = 0; // Value of the biasing function at the current position

    // Energy difference between lattice 1 and lattice 2
    double energy_diff = energy_diff_one_to_two(lattice1_positions, lattice2_positions, params);

    // Adds contributions from all the Gaussians
    for (long gauss_no = 0; gauss_no < params->no_gaussians; gauss_no++) {
        bias_func_val += gaussian(energy_diff, params->gauss_width, params->gauss_positions[gauss_no], params->gauss_heights[gauss_no]);
    }

    double return_val = params->temp * bias_func_val + get_unbiased_energy(current_positions, params);

    if (current_lattice == 1) return_val += params->energy_shift;

    return return_val;
}


// NOT COMPLETED AT ALL
/* Fills an array of forces with the forces on particles at given positions in a certain lattice */
void get_unbiased_forces(double *forces, double *positions, Parameters *params) {
    fill_potential_deriv(forces, positions);
}


void get_biased_forces(int current_lattice, double *bias_forces, double *lattice1_positions, double *lattice2_positions,
                       Parameters *params) {

    fill_double_arr_zeros(bias_forces, params->no_dimensions); // Ensures the biased forces array is empty

    // Finds forces on each lattice
    double *lattice1_forces = malloc(sizeof(double) * params->no_dimensions);
    get_unbiased_forces(lattice1_forces, lattice1_positions, params);

    double *lattice2_forces = malloc(sizeof(double) * params->no_dimensions);
    get_unbiased_forces(lattice2_forces, lattice2_positions, params);

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
    int other_lattice = (params->current_lattice + 1) % 2;

    double current_biased_energy = get_biased_energy(params->current_lattice, lattice1_positions, lattice2_positions,
                                                     params);
    double other_biased_energy   = get_biased_energy(other_lattice, lattice1_positions, lattice2_positions, params);

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

    double *bias_forces = malloc(sizeof(double) * params->no_dimensions); // Array for the bias forces on a particular lattice


    // Decreases f until its at its minimum. At that point the Gaussians are sufficient for our biasing
    while (f > min_f) {
        int flat_histogram = 0; // Indicates if our histogram is flat (0 being not flat, 1 being flat)

        long hist_attempt = 0;
        long tot_timesteps = 0;
        long switch_no = 0;
        long no_left = 0;

        FILE *everything_file = fopen("everything_output.txt", "w");

        // Adds Gaussians until our histogram becomes sufficiently flat
        while (flat_histogram == 0) {
            hist_attempt++;
            if (hist_attempt % 1000 == 0) {
                printf("%ld\n", hist_attempt);

                double dist1 = sqrt(lattice1_positions[0] * lattice1_positions[0] + lattice1_positions[1] * lattice1_positions[1]);
                printf("Distance left from origin = %lf\n", dist1);
            }

            if (hist_attempt % 10000 == 0) {
                printf("Switch_no = %ld, Tot_timesteps = %ld\n", switch_no, tot_timesteps);
                printf("No_left = %ld\n", no_left);
                FILE *bins_file = fopen("bins_output.txt", "w");
                for (long bin_no = 0; bin_no < params->no_bins; bin_no++){
                    fprintf(bins_file, "%ld %ld\n", bin_no, energy_diff_histogram[bin_no]);
                }
                fclose(bins_file);

                FILE *gauss_file = fopen("gauss_output.txt", "w");
                for (long gauss_no = 0; gauss_no < params->no_gaussians; gauss_no++) {
                    fprintf(gauss_file, "%lf\n", params->gauss_positions[gauss_no]);
                }
                fclose(gauss_file);

                fclose(everything_file);

                free(lattice1_positions);
                free(lattice2_positions);
                free(bias_forces);
                free(energy_diff_histogram);
                free_parameters(params);
                exit(0);
            }


            // Performs a set number of BAOAB steps
            for (long stepno = 0; stepno < params->steps_between_gaussians; stepno++) {
                fprintf(everything_file, "%lf %lf %lf %lf %lf\n", lattice1_positions[0], lattice1_positions[1],
                        lattice2_positions[0], lattice2_positions[1], energy_diff_one_to_two(lattice1_positions, lattice2_positions, params));


                // Adds to the histogram with the current energy difference
                add_to_histogram(energy_diff, energy_diff_histogram, params->no_bins, params->minimum_energy_diff,
                                 params->maximum_energy_diff, params->bin_width);

                tot_timesteps++;
                if (params->current_lattice == 0) no_left++;
                // Attempts a lattice switch with set probability
                if (genrand_real1() < params->prob_switch_attempt) {
                    int temp = params->current_lattice;
                    lattice_switch(lattice1_positions, lattice2_positions, params);
                    if (temp != params->current_lattice) switch_no++;
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
            params->gauss_positions[params->no_gaussians] = energy_diff_one_to_two(lattice1_positions, lattice2_positions, params);
            params->gauss_heights[params->no_gaussians] = params->temp * log(f);
            params->no_gaussians++;

            // If the number of Gaussians now exceeds the length of the height and position arrays, extend their length
            if (params->no_gaussians >= params->gauss_arrs_length) {
                long new_arrs_length = params->gauss_arrs_length + params->init_gauss_arrs_length;
                params->gauss_positions = extend_arr(params->gauss_positions, params->gauss_arrs_length, new_arrs_length);
                params->gauss_heights   = extend_arr(params->gauss_heights,   params->gauss_arrs_length, new_arrs_length);
                params->gauss_arrs_length += params->init_gauss_arrs_length;
            }


            // Determines if the histogram is now flat or not
            flat_histogram = is_histogram_flat(energy_diff_histogram, params->no_bins, params->flatness_threshold);
        }
        printf("Got a flat one\n");
        f = sqrt(f); // Decrease f
        fill_long_arr_zeros(energy_diff_histogram, params->no_bins); // Reset the histogram
    }

    free(bias_forces);
    free(energy_diff_histogram);
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
    free(initial_energies);

    implement_biasing(lattice1_positions, lattice2_positions, &params);


    free(lattice1_positions);
    free(lattice2_positions);
    free_parameters(&params);

    return 0;
}