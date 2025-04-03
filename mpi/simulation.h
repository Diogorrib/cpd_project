#ifndef SIMULATION_H
#define SIMULATION_H

#include "constant.h"
#include "globals.h"
#include "particle.h"
#include "utils.h"
#include "comm_utils.h"

void initial_particle_distribution(particle_t *par, cell_t *cells);
void compute_center_of_mass(particle_t *par, cell_t *cells, MPI_Request *r);
void compute_forces(particle_t *par, cell_t *cells, MPI_Request *r);
void compute_forces_maximize(particle_t *par, cell_t *cells, MPI_Request *r);
void compute_new_positions(particle_t **par, cell_t *cells);
long check_collisions(particle_t *par, cell_t *cells);

#endif // SIMULATION_H
