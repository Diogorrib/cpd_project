#ifndef PARTICLE_H
#define PARTICLE_H
#include <mpi.h>


extern MPI_Datatype cell_type;

typedef struct {
    double x, y, vx, vy, m;
    double fx, fy;
    long long cell_idx;
    char is_particle_0;
} particle_t;

typedef struct {
    double x, y, m;
    long long n_part;
    long long *part_idx;
} cell_t;

typedef struct {
    double x, y, m;
} center_of_mass_t;

long long init_particles(long seed, double side, long ncside, long long n_part, long block_low, long block_high, particle_t **first, particle_t **par);

#endif // PARTICLE_H
