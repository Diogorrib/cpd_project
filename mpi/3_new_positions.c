#include "simulation.h"
#include <stdio.h>

void reset_cell_n_parts(cell_t *cells)
{
    for (long long i = 0; i < n_local_cells; i++) { // Last row is ignored by n_local_cells
        long long cell_idx = i + ncside; // Skip the first row as it is computed by another process
        cells[cell_idx].n_part = 0;
    }
}

/**
 * @brief Check if the cell index is from the previous process space
 * 
 * @return true if the cell index is from the previous process space,
 * false if it is from the next process space
 */
int is_prev_process(long long cell_idx)
{
    if (rank == 0 && cell_idx >= (block_size+2) * ncside) return 1; // a lot bigger so it is the last process (previous one with wrap around)
    if (rank == process_count - 1 && cell_idx < 0) return 0; // a lot smaller so it is the first process (next one with wrap around)
    return cell_idx < ncside;
}

/**
 * @brief Ditribute local particles to the process cells and get particles to send to adjacent processes
 */
void distribute_local_parts_and_save_to_exchange(particle_t *par, cell_t *cells, particle_t **prev, particle_t **next, long long *n_prev, long long *n_next)
{
    for (long long i = 0; i < n_part; i++) {
        particle_t *p = &par[i];
        if (p->m == 0) continue;

        long long cell_idx = get_local_cell_idx(p);

        if(cell_in_process_space(cell_idx)) {
            append_particle_to_cell(i, cell_idx, par, cells, 5);
        } else {
            if (is_prev_process(cell_idx)) {
                append_particle_to_array(*n_prev, p, prev, 6);
                (*n_prev)++;
            } else {
                append_particle_to_array(*n_next, p, next, 7);
                (*n_next)++;
            }
            p->m = 0; // mark particle as outside the process space
        }
    }
}

void distribute_received_parts(particle_t *par, cell_t *cells, long long n_part_old)
{
    for (long long i = n_part_old; i < n_part; i++) {
        particle_t *p = &par[i];
        if (p->m == 0) continue;

        long long cell_idx = get_local_cell_idx(p);

        /* if (!cell_in_process_space(cell_idx)) {
            fprintf(stdout, "Rank %d, Error: particle in cell %lld (local=%lld) not in process space; range: %ld - %ld, ncside=%ld\n", rank, get_global_cell_idx(p), cell_idx, block_low*ncside, (block_low+block_size)*ncside-1, ncside);
            exit(EXIT_FAILURE);
        } */

        append_particle_to_cell(i, cell_idx, par, cells, 11);
    }
}

void receive_chuncks_in_loop(particle_t *tmp_prev, particle_t *tmp_next, MPI_Request *prev_req, MPI_Request *next_req, particle_t **local_parts)
{
    int next_count = 0, prev_count = 0;
    int to_recv_next = 1, to_recv_prev = 1;
    int count = 1; // message 0 already being received
    particle_t *tmp_prev_old = NULL;
    particle_t *tmp_next_old = NULL;
    while (to_recv_next || to_recv_prev) {

        // get count to save received particles and start receiving next chunk async while saving current chunk
        if (to_recv_next) {
            next_count = wait_and_get_count(next_req);
            tmp_next_old = tmp_next; // save current chunk reference
            if (next_count >= CHUNK_SIZE) {
                tmp_next = (particle_t *)allocate_memory(CHUNK_SIZE, sizeof(particle_t), 9); // allocate next chunk
                async_recv_part_in_chunks(tmp_next, next_rank, BASE_TAG_1 + 2*count, next_req);
            }
        }

        // get count to save received particles and start receiving next chunk async while saving current chunk
        if (to_recv_prev) {
            prev_count = wait_and_get_count(prev_req);
            tmp_prev_old = tmp_prev; // save current chunk reference
            if (prev_count >= CHUNK_SIZE) {
                tmp_prev = (particle_t *)allocate_memory(CHUNK_SIZE, sizeof(particle_t), 10); // allocate next chunk
                async_recv_part_in_chunks(tmp_prev, prev_rank, BASE_TAG_2 + 2*count, prev_req);
            }
        }

        // save received particles to local array and free tmp array
        if (to_recv_next) convert_to_local_array(tmp_next_old, next_count, local_parts);
        if (to_recv_prev) convert_to_local_array(tmp_prev_old, prev_count, local_parts);

        // last message is smaller than CHUNK_SIZE, others are CHUNK_SIZE
        if (next_count < CHUNK_SIZE) to_recv_next = 0;
        if (prev_count < CHUNK_SIZE) to_recv_prev = 0;
        count++;
    }
}

int exchange_particles(particle_t *prev, particle_t *next, long long n_prev, long long n_next, MPI_Request **requests,
    particle_t *tmp_prev, particle_t *tmp_next, MPI_Request *prev_req, MPI_Request *next_req, particle_t **local_parts)
{
    int n_messages_prev = n_prev / CHUNK_SIZE + 1;
    int n_messages_next = n_next / CHUNK_SIZE + 1;
    int n_messages = n_messages_prev + n_messages_next;

    *requests = (MPI_Request *)allocate_memory(n_messages, sizeof(MPI_Request), 8);
    async_send_part_in_chunks(prev, next, n_prev, n_next, *requests);

    receive_chuncks_in_loop(tmp_prev, tmp_next, prev_req, next_req, local_parts);
    return n_messages;
}

/**
 * @brief Re-distributes particles into cells based on their positions,
 * excluding particles that have collided or outside process space (mass = 0)
 *
 * @param par array of local particles (that will be updated with exchanged particles)
 * @param cells array of local cells accounting for 2 adjacents rows (that will be updated according to the new positions)
 */
void particle_redistribution(particle_t **par, cell_t *cells)
{
    MPI_Request prev_req, next_req;
    MPI_Request *requests;
    particle_t *parts_to_prev = NULL;
    particle_t *parts_to_next = NULL;
    long long n_parts_to_prev = 0;
    long long n_parts_to_next = 0;
    int n_messages;

    particle_t *tmp_prev = (particle_t *)allocate_memory(CHUNK_SIZE, sizeof(particle_t), 3);
    particle_t *tmp_next = (particle_t *)allocate_memory(CHUNK_SIZE, sizeof(particle_t), 4);
    async_recv_part_in_chunks(tmp_next, next_rank, BASE_TAG_1, &next_req);
    async_recv_part_in_chunks(tmp_prev, prev_rank, BASE_TAG_2, &prev_req);

    // reset cells
    reset_cell_n_parts(cells);

    // distribute per cell the particles that are kept in process space
    distribute_local_parts_and_save_to_exchange(*par, cells, &parts_to_prev, &parts_to_next, &n_parts_to_prev, &n_parts_to_next);

    // save current number to distribute received particles
    long long n_part_old = n_part;

    // start sending & wait to receive all particles
    n_messages = exchange_particles(parts_to_prev, parts_to_next, n_parts_to_prev, n_parts_to_next, &requests,
        tmp_prev, tmp_next, &prev_req, &next_req, par);

    // distribute per cell the received particles
    distribute_received_parts(*par, cells, n_part_old);

    wait_for_send_parts(requests, n_messages);
    free(requests);
    free(parts_to_prev);
    free(parts_to_next);
}

void initial_particle_distribution(particle_t *par, cell_t *cells)
{
    // init cells
    reset_cell_n_parts(cells);

    for (long long i = 0; i < n_part; i++) {
        particle_t *p = &par[i];
        if (p->m == 0) continue;

        // here particles should be in process space
        long long cell_idx = get_local_cell_idx(p);
        append_particle_to_cell(i, cell_idx, par, cells, 2);
    }
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
 * check if new positions are outside the space and update the cells.
 * NOTE: Need particle_t** to reallocate the array after receiving new particles
 *
 * @param par array of particles
 * @param cells array of cells
 */
void compute_new_positions(particle_t **par, cell_t *cells)
{
    for (long long i = 0; i < n_part; i++) {
        particle_t *p = &(*par)[i];
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
    particle_redistribution(par, cells);
}
