#include <stdio.h>
#include "comm_utils.h"
#include "utils.h"
#include "globals.h"

#define TAG_CM_1 0 // Tag to send center of mass to previous (the previous receives from next); i.e. send first row and receive ghost row at the end
#define TAG_CM_2 1 // Tag to send center of mass to next (the next receives from previous); i.e. send last row and receive ghost row at the beginning

MPI_Datatype cell_type;
MPI_Datatype part_type;
MPI_Datatype row_type;

/**
 * @brief Receive the center of mass from the adjacent processes asynchronously
 */
void async_recv_center_of_mass(cell_t *cells, MPI_Request *requests)
{
    // Pointers to boundary / adjacent rows
    cell_t *extra_top_row = &cells[0];
    cell_t *extra_bottom_row = &cells[ncside * (block_size + 1)];

    // Receive last (ghost) row from next rank
    MPI_Irecv(extra_bottom_row, 1, row_type, next_rank, TAG_CM_1, MPI_COMM_WORLD, &requests[0]);

    // Receive first (ghost) row from previous rank
    MPI_Irecv(extra_top_row, 1, row_type, prev_rank, TAG_CM_2, MPI_COMM_WORLD, &requests[1]);
}

/**
 * @brief Send the center of mass to the adjacent processes asynchronously
 */
void async_send_center_of_mass(cell_t *cells, MPI_Request *requests)
{
    // Pointers to working rows
    cell_t *first_row = &cells[ncside];
    cell_t *last_row = &cells[ncside * block_size];

    // Send first row to prev rank
    MPI_Isend(first_row, 1, row_type, prev_rank, TAG_CM_1, MPI_COMM_WORLD, &requests[2]);

    // Send last row to next rank
    MPI_Isend(last_row, 1, row_type, next_rank, TAG_CM_2, MPI_COMM_WORLD, &requests[3]);
}

/**
 * @brief Wait for the center of mass to be received and sent
 */
void wait_for_center_of_mass(MPI_Request *requests)//, cell_t *cells)
{
    //MPI_Status status[4];
    MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);//status);

    // argument "cells" needed for debugging
    // debug_cms(status, cells);
}

void async_send_part_in_chunks(particle_t *prev, particle_t *next, long long n_prev, long long n_next, MPI_Request *requests)
{
    int n_messages_prev = n_prev / CHUNK_SIZE + 1;
    int n_messages_next = n_next / CHUNK_SIZE + 1;

    // Send particles to previous rank
    for (int i = 0; i < n_messages_prev; i++) {
        long long size = (i == n_messages_prev - 1) ? n_prev % CHUNK_SIZE : CHUNK_SIZE; // last iteration is smaller
        MPI_Isend(&prev[i * CHUNK_SIZE], size, part_type, prev_rank, BASE_TAG_1 + 2*i, MPI_COMM_WORLD, &requests[i]);
    }

    // Send particles to next rank
    for (int i = 0; i < n_messages_next; i++) {
        long long size = (i == n_messages_next - 1) ? n_next % CHUNK_SIZE : CHUNK_SIZE; // last iteration is smaller
        MPI_Isend(&next[i * CHUNK_SIZE], size, part_type, next_rank, BASE_TAG_2 + 2*i, MPI_COMM_WORLD, &requests[n_messages_prev + i]);
    }
}

void async_recv_part_in_chunks(particle_t *tmp, int rank, int tag, MPI_Request *request)
{
    MPI_Irecv(tmp, CHUNK_SIZE, part_type, rank, tag, MPI_COMM_WORLD, request);
}

int wait_and_get_count(MPI_Request *request)
{
    MPI_Status status;
    int count;

    MPI_Wait(request, &status);

    MPI_Get_count(&status, part_type, &count);

    return count;
}

void wait_for_send_parts(MPI_Request *requests, int count)
{
    MPI_Waitall(count, requests, MPI_STATUSES_IGNORE);
}

void convert_to_local_array(particle_t *tmp, int count, particle_t **local_par)
{
    for (int i = 0; i < count; i++) {
        if (tmp[i].is_particle_0) { // particle_0 received from adjacent process
            //fprintf(stdout, "Rank %d - PARTICLE_0 MOVED HELP!!!", rank);
            particle_0_idx = n_part;
        }

        append_particle_to_array(n_part, &tmp[i], local_par, 12);
        n_part++;
    }
    free(tmp);
}

/**
 * @brief Create MPI data types for the center of mass to be used in communications
 */
void create_mpi_types_for_cms()
{
    int block_lengths[3] = {1, 1, 1};  // Each field is 1 value
    MPI_Aint displacements[3];         // Memory offsets
    MPI_Datatype types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};  // Data types

    displacements[0] = offsetof(cell_t, x);
    displacements[1] = offsetof(cell_t, y);
    displacements[2] = offsetof(cell_t, m);

    // Create the structured datatype
    MPI_Type_create_struct(3, block_lengths, displacements, types, &cell_type);
    MPI_Type_commit(&cell_type);

    // Create a row type to send a full row of cells
    MPI_Type_contiguous(ncside, cell_type, &row_type);
    MPI_Type_commit(&row_type);
}

void create_mpi_particle_type()
{
    int block_lengths[6] = {1, 1, 1, 1, 1, 1};  // Each field is 1 value
    MPI_Aint displacements[6];         // Memory offsets
    MPI_Datatype types[6] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_CHAR};  // Data types
    
    displacements[0] = offsetof(particle_t, x);
    displacements[1] = offsetof(particle_t, y);
    displacements[2] = offsetof(particle_t, vx);
    displacements[3] = offsetof(particle_t, vy);
    displacements[4] = offsetof(particle_t, m);
    displacements[5] = offsetof(particle_t, is_particle_0);

    // Create the structured datatype
    MPI_Type_create_struct(6, block_lengths, displacements, types, &part_type);
    MPI_Type_commit(&part_type);
}

void debug_cms(MPI_Status *status, cell_t *cells)
{
    MPI_Barrier(MPI_COMM_WORLD);
    int count1, count2;
    MPI_Get_count(&status[0], MPI_INT, &count1);
    MPI_Get_count(&status[1], MPI_INT, &count2);
    printf("Rank %d received %d elements from rank %d (tag %d)\n", rank, count1, status[0].MPI_SOURCE, status[0].MPI_TAG);
    printf("Rank %d received %d elements from rank %d (tag %d)\n", rank, count2, status[1].MPI_SOURCE, status[1].MPI_TAG);

    cell_t *extra_top_row = &cells[0];
    cell_t *extra_bottom_row = &cells[ncside * (block_size + 1)];
    cell_t *first_row = &cells[ncside];
    cell_t *last_row = &cells[ncside * block_size];

    for (long long i = 0; i < ncside; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
        fprintf(stdout, "Sent - rank %d to rank %d, cell %lld: x = %.6f, y = %.6f, m = %.6f\n", rank, prev_rank, i, first_row[i].x, first_row[i].y, first_row[i].m);
        MPI_Barrier(MPI_COMM_WORLD);
        fprintf(stdout, "Recv - rank %d from rank %d, cell %lld: x = %.6f, y = %.6f, m = %.6f\n", rank, prev_rank, i, extra_top_row[i].x, extra_top_row[i].y, extra_top_row[i].m);
        MPI_Barrier(MPI_COMM_WORLD);
        fprintf(stdout, "Sent - rank %d to rank %d, cell %lld: x = %.6f, y = %.6f, m = %.6f\n", rank, next_rank, i + ncside * (block_size + 1), last_row[i].x, last_row[i].y, last_row[i].m);
        MPI_Barrier(MPI_COMM_WORLD);
        fprintf(stdout, "Recv - rank %d from rank %d, cell %lld: x = %.6f, y = %.6f, m = %.6f\n", rank, next_rank, i + ncside * (block_size + 1), extra_bottom_row[i].x, extra_bottom_row[i].y, extra_bottom_row[i].m);
        MPI_Barrier(MPI_COMM_WORLD);
    }
}
