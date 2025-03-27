#include "simulation.h"

/**
 * @brief Compute the center of mass for each cell
 * 
 * @param ncside number of cells on each side
 * @param block_size number of cells in the block
 * @param n_part number of particles
 * @param par array of particles
 * @param cells array of cells
 */
void compute_center_of_mass(long ncside, long long block_size, long long n_part, particle_t *par, cell_t *cells)
{
    long long n_cells_process = ncside * block_size;
    (void)n_part;

    for (long long i = 0; i < n_cells_process; i++) {
        cell_t *cell = &cells[i+ncside]; // Skip the first row as it is sent from another process
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
    //TODO: Send the center of mass to the adjacent processes
}
