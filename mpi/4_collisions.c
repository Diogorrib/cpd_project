#include "simulation.h"

long check_collisions_for_part(cell_t *cell, long long j, particle_t *par)
{
    long collisions = 0;
    particle_t *p1 = &par[cell->part_idx[j]];
    for (long long i = j + 1; i < cell->n_part; i++) {
        particle_t *p2 = &par[cell->part_idx[i]];

        double dx = p1->x - p2->x;
        double dy = p1->y - p2->y;
        double dist_sq = dx * dx + dy * dy;
        if (dist_sq <= EPSILON2) {
            if (p1->m != 0 && p2->m != 0) {
                collisions++;
            }
            p1->m = 0;
            p2->m = 0;
        }
    }
    return collisions;
}

/**
 * @brief Check collisions for each particle in each cell,
 * mark particles as collided by setting mass to 0
 *
 * @param par array of particles
 * @param cells array of cells
 * @return long number of collisions
 */
long check_collisions(particle_t *par, cell_t *cells, particle_t *parts_to_prev, particle_t *parts_to_next,
    MPI_Request *requests, int n_messages)
{
    long collisions = 0;

    #pragma omp parallel for schedule(dynamic, 1) reduction(+:collisions)
    for (long long i = 0; i < n_local_cells; i++) { // Last row is ignored by n_local_cells
        long long cell_idx = i + ncside; // Skip the first row as it is computed by another process
        cell_t *cell = &cells[cell_idx];

        for (long long j = 0; j < cell->n_part; j++) {
            collisions += check_collisions_for_part(cell, j, par);
        }
    }

    // wait before next iteration
    wait_for_send_parts(requests, n_messages);
    free(requests);
    free(parts_to_prev);
    free(parts_to_next);

    return collisions;
}
