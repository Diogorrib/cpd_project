#ifndef PARTICLE_H
#define PARTICLE_H

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

long long init_particles(long userseed, double side, long ncside, long long n_part, long long *part_0_idx, particle_t **par);

#endif // PARTICLE_H
