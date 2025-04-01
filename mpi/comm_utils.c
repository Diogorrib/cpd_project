#include <stdio.h>
#include "comm_utils.h"
#include <mpi.h>


void exchange_boundaries(long ncside, long size, cell_t *cells, int rank, int num_procs, MPI_Datatype row_type) {
    MPI_Status status;


    int prev_rank = (rank - 1 + num_procs) % num_procs;
    int next_rank = (rank + 1) % num_procs;


    printf("rank %d exchanging boundaries with rank %d and %d\n", rank, prev_rank, next_rank);

    // Row pointers
    cell_t *first_row = &cells[ncside];                      // First working row
    cell_t *extra_top_row = &cells[0];                       // Extra top row
    cell_t *last_row = &cells[ncside * size];                // Last working row
    cell_t *extra_bottom_row = &cells[ncside * (size + 1)];  // Extra bottom row

    // Send first row to prev rank, receive last row from next rank
    MPI_Sendrecv(first_row, 1, row_type, prev_rank, 0,
                 extra_bottom_row, 1, row_type, next_rank, 0,
                 MPI_COMM_WORLD, &status);            

    // Send last row to next rank, receive first row from prev rank
    MPI_Sendrecv(last_row, 1, row_type, next_rank, 1,
                 extra_top_row, 1, row_type, prev_rank, 1,
                 MPI_COMM_WORLD, &status);
}