#ifndef COMM_UTILS_H
#define COMM_UTILS_H

#include <stdlib.h>
#include <mpi.h>
#include "particle.h"

extern MPI_Datatype cell_type;
extern MPI_Datatype row_type;

void exchange_boundaries(cell_t *cells);
void create_mpi_row_type();
void create_mpi_cell_type();

#endif // COMM_UTILS_H
