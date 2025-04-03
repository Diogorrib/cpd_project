#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "simulation.h"

// define here the globals
double side;
long ncside;
long long n_part;

double inv_cell_side;
long long n_local_cells;
particle_t *particle_0 = NULL;

int rank, process_count;
long block_low, block_size;
int prev_rank, next_rank;

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
long simulation_step(particle_t **par, cell_t *cells, long long time_step)
{   
    (void)time_step;
    MPI_Request center_of_mass_requests[4];
    //fprintf(stdout, "Rank %d - before CM %lld\n", rank, time_step);
    compute_center_of_mass(*par, cells, center_of_mass_requests);
    //fprintf(stdout, "Rank %d - before forces %lld\n", rank, time_step);
    compute_forces(*par, cells, center_of_mass_requests);
    //fprintf(stdout, "Rank %d - before new position %lld\n", rank, time_step);
    compute_new_positions(par, cells);
    //fprintf(stdout, "Rank %d - before collisions %lld\n", rank, time_step);
    return check_collisions(*par, cells);
}
void print_result(long long collisions, double exec_time)
{
    fprintf(stdout, "%.3f %.3f\n", particle_0->x, particle_0->y);
    fprintf(stdout, "%lld\n", collisions);
    fprintf(stderr, "%.1fs\n", exec_time);
}

int main(int argc, char *argv[])
{
    long seed;            // seed for the random number generator
    long long time_steps; // number of time-steps
    long long total_collisions, collisions = 0;

    MPI_Init (&argc, &argv);

    parse_args(argc, argv, &seed, &side, &ncside, &n_part, &time_steps);

    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &process_count);

    block_low = BLOCK_LOW(rank, process_count, ncside);
    block_size = BLOCK_SIZE(rank, process_count, ncside);
    n_local_cells = ncside * block_size;
    inv_cell_side = ncside / side;

    prev_rank = (rank - 1 + process_count) % process_count;
    next_rank = (rank + 1) % process_count;

    // Define mpi types - need to be after MPI_Init and args
    create_mpi_types_for_cms();
    create_mpi_particle_type();

    particle_t *par;
    cell_t *cells = (cell_t *)allocate_memory(ncside*(block_size+2), sizeof(cell_t), 1); // account for adjacent rows
    n_part = init_particles(seed, side, ncside, n_part, &particle_0, &par);

    double exec_time;
    exec_time = -omp_get_wtime();

    initial_particle_distribution(par, cells);
    for (long long t = 0; t < time_steps; t++) {
        collisions += simulation_step(&par, cells, t);
        //fprintf(stdout, "Rank %d - after step %lld\n", rank, t);
    }

    MPI_Allreduce(&collisions, &total_collisions, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD); // already inplicit barrier
    //MPI_Barrier(MPI_COMM_WORLD);

    exec_time += omp_get_wtime();
    if(particle_0 != NULL){
        print_result(total_collisions, exec_time);
    }
    cleanup_cells(cells);
    free(cells);
    free(par);

    MPI_Finalize();
}
