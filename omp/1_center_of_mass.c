#include "simulation.h"
#include <omp.h>

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

    // Initialize cells in parallel
    #pragma omp parallel for
    for (long long i = 0; i < n_cells; i++) {
        cell_t *cell = &cells[i];
        cell->x = 0;
        cell->y = 0;
        cell->m = 0;
    }

    // Accumulate mass and positions in parallel
    #pragma omp parallel for
    for (long long i = 0; i < n_part; i++) {
        particle_t *p = &par[i];
        if (p->m == 0) continue;

        cell_t *cell = &cells[p->cell_idx];

        #pragma omp atomic
        cell->m += p->m;

        #pragma omp atomic
        cell->x += p->m * p->x;

        #pragma omp atomic
        cell->y += p->m * p->y;
    }

    // Compute the center of mass in parallel
    #pragma omp parallel for
    for (long long i = 0; i < n_cells; i++) {
        cell_t *cell = &cells[i];
        if (cell->m != 0) {
            double inv_mass = 1.0 / cell->m;
            cell->x *= inv_mass;
            cell->y *= inv_mass;
        }
    }
}
