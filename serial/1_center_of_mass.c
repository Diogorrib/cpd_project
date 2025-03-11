#include "simulation.h"

/**
 * @brief Compute the center of mass for each cell
 * 
 * @param ncside number of cells on each side
 * @param n_part number of particles
 * @param par array of particles
 * @param cells array of cells
 */
void compute_center_of_mass(long ncside, long long n_part, particle_t *par, cell_t *cells)
{
    long long n_cells = ncside * ncside;
    (void)n_part;

    for (long long i = 0; i < n_cells; i++) {
        cell_t *cell = &cells[i];
        cell->x = 0;
        cell->y = 0;
        cell->m = 0;

        for (long long j = 0; j < cell->n_part; j++) {
            particle_t *p = &par[cell->part_idx[j]];
            cell->m += p->m;
            cell->x += p->m * p->x;
            cell->y += p->m * p->y;
        }

        if (cell->m != 0) {
            double inv_mass = 1.0 / cell->m;
            cell->x *= inv_mass;
            cell->y *= inv_mass;
        }
    }
}
