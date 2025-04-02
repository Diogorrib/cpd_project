#include "simulation.h"
#include <stdio.h>

/**
 * @brief Distributes particles into cells based on their positions,
 * excluding particles that have collided (mass = 0)
 *
 * @param par array of particles
 * @param cells array of cells
 */
void particle_distribution(particle_t *par, cell_t *cells)
{
    particle_t *parts_to_prev = NULL;
    particle_t *parts_to_next = NULL;
    long long n_parts_to_prev = 0;
    long long n_parts_to_next = 0;

    MPI_Request prev_req;
    MPI_Request next_req;

    particle_t *tmp_prev = malloc(CHUNK_SIZE * sizeof(particle_t));
    particle_t *tmp_next = malloc(CHUNK_SIZE * sizeof(particle_t));
    if (!tmp_prev || !tmp_next) {
        fprintf(stderr, "Rank %d, Memory allocation failed (10)", rank);
        exit(EXIT_FAILURE);
    }
    async_recv_part_in_chunks(tmp_next, next_rank, BASE_TAG_1, &next_req);
    async_recv_part_in_chunks(tmp_prev, prev_rank, BASE_TAG_2, &prev_req);

    // init / reset cells
    for (long long i = 0; i < n_local_cells; i++) { // Last row is ignored by n_local_cells
        long long cell_idx = i + ncside; // Skip the first row as it is computed by another process
        cells[cell_idx].n_part = 0;
    }

    for (long long i = 0; i < n_part; i++) {
        particle_t *p = &par[i];
        if (p->m == 0) continue;

        long long cell_idx = get_local_cell_idx(p);

        if(cell_in_process_space(cell_idx)) {
            append_particle_to_cell(i, cell_idx, par, cells);
        } else {
            if (cell_idx < ncside) {
                append_particle_to_array(n_parts_to_prev, p, &parts_to_prev);
                n_parts_to_prev++;
            } else {
                append_particle_to_array(n_parts_to_next, p, &parts_to_next);
                n_parts_to_next++;
            }
            p->m = 0; // mark particle as outside the process space
        }
        /* if (i == n_part - 1 && n_parts_to_prev > 0) fprintf(stdout, "Rank %d: %lld particles to prev rank\n", rank, n_parts_to_prev);
        if (i == n_part - 1 && n_parts_to_next > 0) fprintf(stdout, "Rank %d: %lld particles to next rank\n", rank, n_parts_to_next); */
    }

    int n_messages = n_parts_to_prev / CHUNK_SIZE + 1 + n_parts_to_next / CHUNK_SIZE + 1;
    MPI_Request *requests = malloc(n_messages * sizeof(MPI_Request));
    if (!requests) {
        fprintf(stderr, "Rank %d, Memory allocation failed (11)", rank);
        exit(EXIT_FAILURE);
    }
    async_send_part_in_chunks(parts_to_prev, parts_to_next, n_parts_to_prev, n_parts_to_next, requests);

    int next_count = 0, prev_count = 0;
    int to_recv_next = 1, to_recv_prev = 1;
    int count = 0;
    particle_t *tmp_prev_old = NULL;
    particle_t *tmp_next_old = NULL;
    while (to_recv_next || to_recv_prev) {
        if (to_recv_next) {
            next_count = wait_and_get_count(&next_req);
            tmp_next_old = tmp_next;
            if (next_count >= CHUNK_SIZE) {
                tmp_next = malloc(CHUNK_SIZE * sizeof(particle_t));
                if (!tmp_next) {
                    fprintf(stderr, "Rank %d, Memory allocation failed (12)", rank);
                    exit(EXIT_FAILURE);
                }
                async_recv_part_in_chunks(tmp_next, next_rank, BASE_TAG_1 + 2*count, &next_req);
            }
        }

        if (to_recv_prev) { 
            prev_count = wait_and_get_count(&prev_req);
            tmp_prev_old = tmp_prev;
            if (prev_count >= CHUNK_SIZE) {
                tmp_prev = malloc(CHUNK_SIZE * sizeof(particle_t));
                if (!tmp_next) {
                    fprintf(stderr, "Rank %d, Memory allocation failed (13)", rank);
                    exit(EXIT_FAILURE);
                }
                async_recv_part_in_chunks(tmp_prev, prev_rank, BASE_TAG_2 + 2*count, &prev_req);
            }
        }

        if (to_recv_next) convert_to_local_array(tmp_next_old, next_count, &par);
        if (to_recv_prev) convert_to_local_array(tmp_prev_old, prev_count, &par);

        if (next_count < CHUNK_SIZE) to_recv_next = 0;
        if (prev_count < CHUNK_SIZE) to_recv_prev = 0;
        count++;
    }

    wait_for_send_parts(requests, n_messages);
    free(requests);
}

void check_outside_space(particle_t *p)
{
    if (p->x < 0) {
        p->x += side;
    } else if (p->x >= side) {
        p->x -= side;
    }

    if (p->y < 0) {
        p->y += side;
    } else if (p->y >= side) {
        p->y -= side;
    }
}

/**
 * @brief  Compute the new velocity and position for each particle,
 * check if new positions are outside the space and update the cells
 *
 * @param par array of particles
 * @param cells array of cells
 */
void compute_new_positions(particle_t *par, cell_t *cells)
{
    for (long long i = 0; i < n_part; i++) {
        particle_t *p = &par[i];
        if (p->m == 0) continue;

        double inv_mass = 1.0 / p->m;
        double ax = p->fx * inv_mass;
        double ay = p->fy * inv_mass;

        p->x += p->vx * DELTAT + 0.5 * ax * DELTAT * DELTAT;
        p->y += p->vy * DELTAT + 0.5 * ay * DELTAT * DELTAT;
        p->vx += ax * DELTAT;
        p->vy += ay * DELTAT;
        check_outside_space(p);

        // reset force for next iteration
        p->fx = 0;
        p->fy = 0;
    }
    cleanup_cells(cells);
    particle_distribution(par, cells);
}
