#include "simulation.h"
#include <mpi.h>

/**
 * @brief Compute the center of mass for each cell
 *
 * @param par array of particles
 * @param cells array of cells
 */
void compute_center_of_mass(particle_t *par, cell_t *cells, MPI_Request *r)
{
    async_recv_center_of_mass(cells, r);

    for (long long i = 0; i < n_local_cells; i++) { // Last row is ignored by n_local_cells
        long long cell_idx = i + ncside; // Skip the first row as it is computed by another process
        cell_t *cell = &cells[cell_idx];
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

    async_send_center_of_mass(cells, r);
}
