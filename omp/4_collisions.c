#include "simulation.h"

long check_collisions_for_part(cell_t *cell, long long j, particle_t *par)
{
    long collisions = 0;
    particle_t *p1 = &par[cell->part_idx[j]];

    for (long long i = j + 1; i < cell->n_part; i++) {
        particle_t *p2 = &par[cell->part_idx[i]];
        if (p1->m == 0 && p2->m == 0) continue;
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
 * @param ncside number of cells on each side
 * @param par array of particles
 * @param cells array of cells
 * @return long number of collisions
 */
long check_collisions(long ncside, particle_t *par, cell_t *cells)
{
    long collisions = 0;
    long long n_cells = ncside * ncside;

    #pragma omp parallel for reduction(+:collisions) schedule(guided, 1)
    for (long long i = 0; i < n_cells; i++) {
        cell_t *cell = &cells[i];
        for (long long j = 0; j < cell->n_part; j++) {
            collisions += check_collisions_for_part(cell, j, par);
        }
    }
    return collisions;
}
