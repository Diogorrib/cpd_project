#include <stdio.h>
#include "comm_utils.h"
#include "globals.h"

MPI_Datatype cell_type;
MPI_Datatype row_type;

void exchange_boundaries(cell_t *cells)
{
    MPI_Status status;

    int prev_rank = (rank - 1 + process_count) % process_count;
    int next_rank = (rank + 1) % process_count;

    printf("rank %d exchanging boundaries with rank %d and %d\n", rank, prev_rank, next_rank);

    // Row pointers
    cell_t *first_row = &cells[ncside];                      // First working row
    cell_t *extra_top_row = &cells[0];                       // Extra top row
    cell_t *last_row = &cells[ncside * block_size];                // Last working row
    cell_t *extra_bottom_row = &cells[ncside * (block_size + 1)];  // Extra bottom row

    // Send first row to prev rank, receive last row from next rank
    MPI_Sendrecv(first_row, 1, row_type, prev_rank, 0,
                 extra_bottom_row, 1, row_type, next_rank, 0,
                 MPI_COMM_WORLD, &status);

    // Send last row to next rank, receive first row from prev rank
    MPI_Sendrecv(last_row, 1, row_type, next_rank, 1,
                 extra_top_row, 1, row_type, prev_rank, 1,
                 MPI_COMM_WORLD, &status);
}

void create_mpi_row_type() 
{
    MPI_Type_contiguous(ncside, MPI_DOUBLE, &row_type);
    MPI_Type_commit(&row_type);
}

void create_mpi_cell_type() 
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
}
