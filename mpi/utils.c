#include <stdio.h>
#include <mpi.h>
#include "utils.h"

MPI_Datatype cell_type;

void* allocate_memory(size_t n_elements, size_t element_size)
{
    void *array = malloc(n_elements * element_size);
    if (!array) {
        fprintf(stderr, "Memory allocation failed (2)\n");
        exit(EXIT_FAILURE);
    }
    return array;
}

long long get_dynamic_chunk_size(long long n_part)
{
    if (n_part < 100) {
        return 10;
    } else if (n_part < 250) {
        return 25;
    } else if (n_part < 500) {
        return 50;
    } else if (n_part < 1000) {
        return 100;
    } else if (n_part < 2500) {
        return 250;
    } else if (n_part < 5000) {
        return 500;
    } else if (n_part < 100000) {
        return 1000;
    } else {
        return 2500;
    }
}

void append_particle_index(long long idx, long long cell_idx, particle_t *par, cell_t *cells)
{
    cell_t *cell = &cells[cell_idx];
    long long n_part = cell->n_part;
    long long *part_idx = cell->part_idx;

    long long chunk_size = get_dynamic_chunk_size(n_part);
    if (n_part == 0) {
        part_idx = (long long *)malloc(chunk_size * sizeof(long long));
    } else if (n_part % chunk_size == 0) {
        part_idx = (long long *)realloc(part_idx, (n_part + chunk_size) * sizeof(long long));
    }
    if (!part_idx) {
        fprintf(stderr, "Memory allocation failed (3) %lld\n", cell_idx);
        exit(EXIT_FAILURE);
    }

    part_idx[n_part] = idx;
    cell->n_part++;
    cell->part_idx = part_idx;

    par[idx].cell_idx = cell_idx;
}

/* int part_process_space(long ncside, double side, long block_low, long block_high, particle_t *p){
    double inv_cell_side = ncside / side;
    long long cell_x_idx = get_cell_idx(inv_cell_side, ncside, p);
    if(cell_x_idx >= block_low*ncside && cell_x_idx < block_high*ncside){
        return 1;
    }
    return 0;
} */

/**
 * @brief Check if the cell index is in the process space
 */
int cell_process_space(long ncside, long block_size, long long cell_idx){
    long long aux_idx = cell_idx - ncside;
    if(aux_idx >= 0 && aux_idx < block_size*ncside) {
        return 1;
    }
    return 0;
}

/**
 * @brief Get local cell index for a particle
 */
long long get_cell_idx(double inv_cell_side, long ncside, long block_low, particle_t *p){
    long cell_x_idx = (long)(p->x * inv_cell_side);
    long cell_y_idx = (long)(p->y * inv_cell_side);
    long long global_cell_idx = cell_x_idx + cell_y_idx * ncside;
    return global_cell_idx - (block_low-1) * ncside;
}

void cleanup_cells(long ncside, long long blocks_size, cell_t *cells)
{
    for (long long i = 0; i < (ncside+2) * blocks_size; i++) {
        cell_t *cell = &cells[i];
        if (cell->n_part > 0) {
            free(cell->part_idx);
        }
    }
}

void create_mpi_cell_type() {
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
