#include "simulation.h"
#include <mpi.h>

/**
 * @brief Compute the center of mass for each cell
 *
 * @param par array of local particles
 * @param cells array of local cells accounting for 2 adjacents rows
 * @param r array of MPI requests to start receiving and to start sending
 */
void compute_center_of_mass(particle_t *par, cell_t *cells, MPI_Request *r)
{
    #pragma omp parallel for schedule(dynamic, 1)
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

void async_start_recv_all(cell_t *cells, MPI_Request *requests_cm, particle_t **tmp_prev, particle_t **tmp_next,
    MPI_Request *prev_req, MPI_Request *next_req)
{
    async_recv_center_of_mass(cells, requests_cm);

    *tmp_prev = (particle_t *)allocate_memory(CHUNK_SIZE, sizeof(particle_t), 3);
    *tmp_next = (particle_t *)allocate_memory(CHUNK_SIZE, sizeof(particle_t), 4);
    async_recv_part_in_chunks(*tmp_next, next_rank, BASE_TAG_1, next_req);
    async_recv_part_in_chunks(*tmp_prev, prev_rank, BASE_TAG_2, prev_req);
}
