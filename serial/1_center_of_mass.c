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

    for (long long i = 0; i < n_cells; i++) {
        cells[i].x = 0;
        cells[i].y = 0;
        cells[i].m = 0;
    }

    for (long long i = 0; i < n_part; i++) {
        if (par[i].m == 0) {
            continue;
        }
        long long cell_idx = par[i].cell_idx;
        cells[cell_idx].m += par[i].m;
        cells[cell_idx].x += par[i].m * par[i].x;
        cells[cell_idx].y += par[i].m * par[i].y;
    }

    for (long long i = 0; i < n_cells; i++) {
        if (cells[i].m != 0) {
            cells[i].x /= cells[i].m;
            cells[i].y /= cells[i].m;
        }
    }
}
