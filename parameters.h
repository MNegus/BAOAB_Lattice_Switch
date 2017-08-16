/*
 * parameters.h
 *
 * Header file for all functions dealing with reading and storing relevant parameters needed. To use, create an input
 * file containing your desired parameters on separate lines in this order:
 *
 * no_dimensions              = Number of dimensions (i.e. 3 * N for N particles in three spatial dimensions)
 * temp                    = Temperature given in units of kT (k = Boltzmann constant, T = temperature in Kelvin)
 * mass                    = Mass of particles (in kg)
 * timestep                = Timestep size between each BAOAB step (seconds)
 * prob_switch_attempt     = Probability of attempting a lattice switch as opposed to a BAOAB step at any timestep
 * initial_gauss_height    = Height of the first Gaussian to be placed down in biasing
 * min_gauss_height        = The minimum height of a Gaussian when bias
 * gauss_width             = Width of the Gaussians
 * steps_between_gaussians = Number of timesteps between each Gaussian being placed
 * minimum_energy_diff     = Minimum of the energy differences between the two lattices
 * maximum_energy_diff     = Maximum of the energy differences between the two lattices
 * no_bins                 = Number of bins in the histogram for recording the energy differences
 * flatness_threshold      = When every bin has greater than flatness_threshold * mean_of_bins of hits, the histogram is flat
 * start_lattice           = Lattice to start the simulations in (either 0 or 1)
 * tot_timesteps           = Total number of timesteps to run the final simulation for
 *
 * Then, create an instance of the Parameters structure, and call store_parameters, giving it a pointer to your instance,
 * the name of your input file, and a 2 element array of the energies of the lattices at their initial configuration,
 * which should be when they are at their minima.
 *
 */

#ifndef BAOAB_LATTICE_SWITCH_PARAMETERS_H
#define BAOAB_LATTICE_SWITCH_PARAMETERS_H

#include "miscfunctions.h"


/* Structure for storing parameters needed for the simulation */
struct Parameters {
    /* Parameters to be read in from an input file */
    // Physical parameters
    long no_dimensions; // Number of no_dimensions in the problem, i.e. 3N for N particles in 3 no_dimensions
    double temp; // Temperature in units of energy, equal to kT, where k = Boltzmann constant, T = temperature in Kelvin
    double mass; // Mass of each particle
    double timestep; // Size of the timestep to be used in the BAOAB limit method
    double prob_switch_attempt; // The probability of attempting a lattice switch as opposed to a BAOAB step at any timestep

    // Biasing parameters
    double initial_gauss_height; // Height of the first Gaussian to be put down
    double min_gauss_height; // Minimum height for the Gaussians
    double gauss_width; // Width of the Gaussians
    long steps_between_gaussians; // BAOAB steps between each Gaussian being placed
    double minimum_energy_diff; // Minimum of the energy differences between the two lattices
    double maximum_energy_diff; // Maximum between the energy differences between the two lattices
    long no_bins; // Number of bins in histogram for recording energy differences
    double flatness_threshold; // When every bin has greater than flatness_threshold * mean_of_bins of hits, the histogram is "flat"

    // Parameters for final simulation, once biasing has been applied
    int start_lattice; // The lattice to start the dynamics in (either 0 or 1)
    long tot_timesteps; // Total number of timesteps to run simulation for


    /* Data structures and variables derived using the above parameters */
    double bin_width; // Width of each bin in the histogram
    double energy_shift; // Amount to shift lattice 1 by to make its minimum energy equal to lattice 0
    double *normal_dist_arr; // Array of normally distributed numbers, length "no_dimensions", used in BAOAB steps
};
typedef struct Parameters Parameters;


/* Reads required parameters from an input file and stores them in a given instance of Parameters */
void store_parameters(Parameters *params, char *input_filename, double *initial_energies) {
    // Attempts to open the input file
    FILE *input_file = fopen(input_filename, "r");
    if (input_file == NULL) {
        printf("Failed to open input file\n");
        exit(1);
    }

    // Reads the input file and stores them into the given parameters struct
    read_long(input_file,   &params->no_dimensions);
    read_double(input_file, &params->temp);
    read_double(input_file, &params->mass);
    read_double(input_file, &params->timestep);
    read_double(input_file, &params->prob_switch_attempt);

    read_double(input_file, &params->initial_gauss_height);
    read_double(input_file, &params->min_gauss_height);
    read_double(input_file, &params->gauss_width);
    read_long(input_file,   &params->steps_between_gaussians);
    read_double(input_file, &params->minimum_energy_diff);
    read_double(input_file, &params->maximum_energy_diff);
    read_long(input_file,   &params->no_bins);
    read_double(input_file, &params->flatness_threshold);

    read_int(input_file,    &params->start_lattice);
    read_long(input_file,   &params->tot_timesteps);

    fclose(input_file);


    // Calculates derived parameters
    params->bin_width = (params->maximum_energy_diff - params->minimum_energy_diff) / params->no_bins;
    params->energy_shift = initial_energies[0] - initial_energies[1]; // The amount lattice 1 has to be shifted

    params->normal_dist_arr = malloc(sizeof(double) * params->no_dimensions); // Ensure the seed has been set before this
    fill_arr_normal_dist(params->normal_dist_arr, params->no_dimensions); // Fills the array with normally distributed numbers
}


/* Frees all manually allocated memory used by a Parameters instance */
void free_parameters(Parameters *params) {
    free(params->normal_dist_arr);
}


#endif //BAOAB_LATTICE_SWITCH_PARAMETERS_H
