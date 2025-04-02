#ifndef COMM_UTILS_H
#define COMM_UTILS_H

#include <stdlib.h>
#include <mpi.h>
#include "particle.h"

#define CHUNK_SIZE 500 // Size of the chunks to send
#define BASE_TAG_1 2 // Base tag to send particles (it cannot be the same because there could be multiple chunks)
#define BASE_TAG_2 3 // Base tag to send particles (it cannot be the same because there could be multiple chunks)

extern MPI_Datatype cell_type;
extern MPI_Datatype part_type;
extern MPI_Datatype row_type;

void async_recv_center_of_mass(cell_t *cells, MPI_Request *requests);
void async_send_center_of_mass(cell_t *cells, MPI_Request *requests);
void wait_for_center_of_mass(MPI_Request *requests);//, cell_t *cells);

void async_send_part_in_chunks(particle_t *prev, particle_t *next, long long n_prev, long long n_next, MPI_Request *requests);
void async_recv_part_in_chunks(particle_t *tmp, int rank, int tag, MPI_Request *request);
int wait_and_get_count(MPI_Request *request);
void wait_for_send_parts(MPI_Request *requests, int count);
void convert_to_local_array(particle_t *tmp, int count, particle_t **local_par);

void create_mpi_types_for_cms();
void create_mpi_particle_type();

#endif // COMM_UTILS_H
