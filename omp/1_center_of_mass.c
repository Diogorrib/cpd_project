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
    int num_threads = omp_get_max_threads();

    cell_t *thread_cells = (cell_t *)calloc(num_threads * n_cells, sizeof(cell_t));

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        cell_t *local_cells = &thread_cells[tid * n_cells];

        #pragma omp for 
        for (long long i = 0; i < n_cells; i++) {
            cells[i].x = 0;
            cells[i].y = 0;
            cells[i].m = 0;
        }

        // Accumulate mass and weighted positions into thread-local storage
        #pragma omp for 
        for (long long i = 0; i < n_part; i++) {
            particle_t *p = &par[i];
            if (p->m == 0) continue;

            cell_t *cell = &local_cells[p->cell_idx];
            cell->m += p->m;
            cell->x += p->m * p->x;
            cell->y += p->m * p->y;
        }

        #pragma omp for
        for (long long i = 0; i < n_cells; i++) {
            for (int t = 0; t < num_threads; t++) {
                cells[i].m += thread_cells[t * n_cells + i].m;
                cells[i].x += thread_cells[t * n_cells + i].x;
                cells[i].y += thread_cells[t * n_cells + i].y;
            }
        }

        #pragma omp for
        for (long long i = 0; i < n_cells; i++) {
            if (cells[i].m != 0) {
                double inv_mass = 1.0 / cells[i].m;
                cells[i].x *= inv_mass;
                cells[i].y *= inv_mass;
            }
        }
    }

    free(thread_cells);
}
