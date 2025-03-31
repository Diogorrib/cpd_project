#ifndef COMM_UTILS_H
#define COMM_UTILS_H

#include <stdlib.h>
#include <mpi.h>
#include "particle.h"

void exchange_boundaries(long ncside, long size, cell_t *cells, int rank, int num_procs, MPI_Datatype row_type);

#endif