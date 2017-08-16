//
// Created by michael on 16/08/17.
//

#ifndef BAOAB_LATTICE_SWITCH_PARAMETERS_H
#define BAOAB_LATTICE_SWITCH_PARAMETERS_H

struct Parameters {

    /* Parameters to be read in from an input file */
    // Physical parameters
    long dimensions; // Number of dimensions in the problem, i.e. 3N for N particles in 3 dimensions
    double temp; // Temperature in units of energy, equal to kT, where k = Boltzmann constant, T = temperature in Kelvin
    double *masses; // Array of length "dimensions" containing the masses for each particle
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
    double *normal_dist_arr; // Array of normally distributed numbers, length "dimensions", used in BAOAB steps
    double energy_shift; // Amount to shift lattice 1 by to make its minimum energy equal to lattice 0
};
typedef struct Parameters Parameters;

#endif //BAOAB_LATTICE_SWITCH_PARAMETERS_H
