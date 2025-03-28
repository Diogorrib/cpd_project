#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>
#include "simulation.h"

#define BLOCK_LOW(id, p, ncside) ((id) * (ncside) / (p))
#define BLOCK_HIGH(id, p, ncside) (BLOCK_LOW((id) + 1, p, ncside) - 1)
#define BLOCK_SIZE(id, p, ncside) (BLOCK_HIGH(id, p, ncside) - BLOCK_LOW(id, p, ncside) + 1)
#define BLOCK_OWNER(index, p, ncside) (((p) * ((index) + 1) - 1) / (ncside))

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
    *time_steps = atoll(argv[5]);
}
long simulation_step(double side, long ncside, long long block_size, long block_low, long long n_part, particle_t *par, cell_t *cells)
{
    compute_center_of_mass(ncside, block_size, par, cells);
    compute_forces(side, ncside, block_size, par, cells);
    compute_new_positions(side, ncside, block_size, block_low, n_part, par, cells);
    return check_collisions(ncside, block_size, par, cells);
}
void print_result(particle_t *part_0, long long collisions, double exec_time)
{
    fprintf(stdout, "%.3f %.3f\n", part_0->x, part_0->y);
    fprintf(stdout, "%lld\n", collisions);
    fprintf(stderr, "%.1fs\n", exec_time);
}

int main(int argc, char *argv[])
{
    long seed;              // seed for the random number generator
    double side;            // size of the side of the squared space of simulation
    long ncside;            // size of the grid (number of cells on each side)ls
    long long n_part;       // number of particles
    long long time_steps;   // number of time-steps
    long long total_collisions, collisions = 0;

    int rank, processes_count = 0;

    long block_low, block_high, block_size;

    MPI_Init (&argc, &argv);

    parse_args(argc, argv, &seed, &side, &ncside, &n_part, &time_steps);

    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &processes_count);

    block_low = BLOCK_LOW(rank, processes_count, ncside);
    block_high = BLOCK_HIGH(rank, processes_count, ncside);
    block_size = BLOCK_SIZE(rank, processes_count, ncside);

    /* fprintf(stdout, "Rank: %d, Block_low: %ld, Block_high: %ld, Block_size: %ld\n",
        rank, block_low, block_high, block_size); */

    particle_t *par;
    particle_t *particle_0 = NULL;
    cell_t *cells = (cell_t *)allocate_memory(ncside*(block_size+2), sizeof(cell_t)); // account for adjacent rows
    n_part = init_particles(seed, side, ncside, n_part, block_low, block_high, &particle_0, &par);

    double exec_time;
    exec_time = -omp_get_wtime();

    particle_distribution(side, ncside, block_size, block_low, n_part, par, cells);
    for (long long t = 0; t < time_steps; t++) {
        collisions += simulation_step(side, ncside, block_size, block_low, n_part, par, cells);
    }

    MPI_Allreduce(&collisions, &total_collisions, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    exec_time += omp_get_wtime();
    if(particle_0 != NULL){
        print_result(particle_0, total_collisions, exec_time); 
    }
    cleanup_cells(ncside, block_size, cells);
    free(cells);
    free(par);

    MPI_Finalize();
}
