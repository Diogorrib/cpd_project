#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "simulation.h"

void parse_args(int argc, char *argv[], long *seed, double *side, long *ncside, long long *n_part, long long *time_steps)
{
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <seed> <side> <ncside> <n_part> <time_steps>\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    *seed = atol(argv[1]);
    *side = atof(argv[2]);
    *ncside = atol(argv[3]);
    *n_part = atoll(argv[4]);
    *time_steps = atol(argv[5]);
}
long simulation_step(double side, long ncside, long long n_part, particle_t *par, cell_t *cells)
{
    compute_center_of_mass(ncside, n_part, par, cells);
    compute_forces(side, ncside, par, cells);
    compute_new_positions(side, ncside, n_part, par, cells);
    return check_collisions(ncside, par, cells);
}
void print_result(particle_t *par, long long collisions)
{
    fprintf(stdout, "%.3f %.3f\n", par[0].x, par[0].y);
    fprintf(stdout, "%lld\n", collisions);
}

int main(int argc, char *argv[])
{
    long seed;              // seed for the random number generator
    double side;            // size of the side of the squared space of simulation
    long ncside;            // size of the grid (number of cells on each side)
    long long n_part;       // number of particles
    long long time_steps;   // number of time-steps
    long long collisions = 0;

    parse_args(argc, argv, &seed, &side, &ncside, &n_part, &time_steps);

    particle_t *par = (particle_t *)allocate_memory(n_part, sizeof(particle_t));
    cell_t *cells = (cell_t *)allocate_memory(ncside * ncside, sizeof(cell_t));
    init_particles(seed, side, ncside, n_part, par);

    double exec_time;
    exec_time = -omp_get_wtime();

    particle_distribution(side, ncside, n_part, par, cells);
    for (long long t = 0; t < time_steps; t++) {
        collisions += simulation_step(side, ncside, n_part, par, cells);
    }

    exec_time += omp_get_wtime();
    fprintf(stderr, "%.1fs\n", exec_time);

    print_result(par, collisions);
    cleanup_cells(ncside, cells);
    free(cells);
    free(par);
}
