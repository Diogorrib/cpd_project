#ifndef SIMULATION_H
#define SIMULATION_H

#include "constant.h"
#include "particle.h"
#include "utils.h"

void particle_distribution(double side, long ncside, long block_size, long block_low, long long n_part, particle_t *par, cell_t *cells);
void compute_center_of_mass(long ncside, long block_size, particle_t *par, cell_t *cells);
void compute_forces(double side, long ncside, long block_size, particle_t *par, cell_t *cells);
void compute_new_positions(double side, long ncside, long block_size, long block_low, long long n_part, particle_t *par, cell_t *cells);
long check_collisions(long ncside, long block_size, particle_t *par, cell_t *cells);

#endif // SIMULATION_H
