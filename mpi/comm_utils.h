#ifndef COMM_UTILS_H
#define COMM_UTILS_H

#include <stdlib.h>
#include <mpi.h>
#include "particle.h"

extern MPI_Datatype cell_type;
extern MPI_Datatype row_type;

void async_recv_center_of_mass(cell_t *cells, MPI_Request *requests);
void async_send_center_of_mass(cell_t *cells, MPI_Request *requests);
void wait_for_center_of_mass(MPI_Request *requests);//, cell_t *cells);
void create_mpi_types_for_cms();

#endif // COMM_UTILS_H
