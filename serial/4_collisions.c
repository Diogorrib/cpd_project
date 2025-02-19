#include "simulation.h"

long check_collisions_for_part(long long cell_idx, long long j, particle_t *par, cell_t *cells)
{
    long collisions = 0;
    long long p1_idx = cells[cell_idx].part_idx[j];
    for (long long i = j + 1; i < cells[cell_idx].n_part; i++) {
        long long p2_idx = cells[cell_idx].part_idx[i];
        double dx = par[p1_idx].x - par[p2_idx].x;
        double dy = par[p1_idx].y - par[p2_idx].y;
        double dist_sq = dx * dx + dy * dy;
        if (dist_sq <= EPSILON2) {
            par[p1_idx].m = 0;
            par[p2_idx].m = 0;
            collisions++;
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
    for (long long i = 0; i < n_cells; i++) {
        for (long long j = 0; j < cells[i].n_part; j++) {
            long long p_idx = cells[i].part_idx[j];
            if (par[p_idx].m != 0) {
                collisions += check_collisions_for_part(i, j, par, cells);
            }
        }
    }
    return collisions;
}
