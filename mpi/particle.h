#ifndef PARTICLE_H
#define PARTICLE_H

typedef struct {
    double x, y, vx, vy, m;
    double fx, fy;
    long long cell_idx;
} particle_t;

typedef struct {
    double x, y, m;
    long long n_part;
    long long *part_idx;
} cell_t;

typedef struct {
    double x, y, m;
} center_of_mass_t;

void init_particles(long seed, double side, long ncside, long long n_part, particle_t *par);

#endif // PARTICLE_H
