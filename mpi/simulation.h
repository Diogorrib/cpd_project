#ifndef SIMULATION_H
#define SIMULATION_H

#include "constant.h"
#include "globals.h"
#include "particle.h"
#include "utils.h"
#include "comm_utils.h"

void initial_particle_distribution(particle_t *par, cell_t *cells);

void async_start_recv_all(cell_t *cells, MPI_Request *requests_cm, particle_t **tmp_prev, particle_t **tmp_next,
    MPI_Request *prev_req, MPI_Request *next_req);

void compute_center_of_mass(particle_t *par, cell_t *cells, MPI_Request *r);

void compute_forces(particle_t *par, cell_t *cells, MPI_Request *r);

int compute_new_positions(particle_t **par, cell_t *cells, particle_t *tmp_prev, particle_t *tmp_next,
    MPI_Request *prev_req, MPI_Request *next_req, particle_t **parts_to_prev, particle_t **parts_to_next,
    MPI_Request **requests);

long check_collisions(particle_t *par, cell_t *cells, particle_t *parts_to_prev, particle_t *parts_to_next,
    MPI_Request *requests, int n_messages);

#endif // SIMULATION_H
