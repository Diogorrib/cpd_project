#include <stdio.h>
#include "utils.h"
#include "globals.h"
#include <mpi.h>

void* allocate_memory(size_t n_elements, size_t element_size)
{
    void *array = calloc(n_elements, element_size);
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

void append_particle_to_array(long long idx, particle_t *p, particle_t **p_array)
{
    long long chunk_size = get_dynamic_chunk_size(idx);

    if (idx == 0) {
        *p_array = (particle_t *)malloc(chunk_size * sizeof(particle_t));
    } else if (idx % chunk_size == 0) {
        *p_array = (particle_t *)realloc(*p_array, (idx + chunk_size) * sizeof(particle_t));
    }
    if (!*p_array) {
        fprintf(stderr, "Memory allocation failed (1)\n");
        exit(EXIT_FAILURE);
    }

    // useful for init & receiving
    p->fx = 0;
    p->fy = 0;

    // assign particle to the array
    (*p_array)[idx] = *p;
}

void append_particle_to_cell(long long idx, long long cell_idx, particle_t *par, cell_t *cells)
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

int cell_in_process_space(long long cell_idx)
{
    long long aux_idx = cell_idx - ncside;
    if (aux_idx >= 0 && aux_idx < block_size * ncside) {
        return 1;
    }
    return 0;
}

long long get_global_cell_idx(particle_t *p)
{
    long cell_x_idx = (long)(p->x * inv_cell_side);
    long cell_y_idx = (long)(p->y * inv_cell_side);
    return cell_x_idx + cell_y_idx * ncside;
}

long long get_local_cell_idx(particle_t *p)
{
    long long global_idx = get_global_cell_idx(p);
    return global_idx - (block_low - 1) * ncside;
}

void cleanup_cells(cell_t *cells)
{
    for (long long i = 0; i < n_local_cells; i++) { // Last row is ignored by n_local_cells
        long long cell_idx = i + ncside; // Skip the first row as it is computed by another process
        cell_t *cell = &cells[cell_idx];
        if (cell->n_part > 0) {
            free(cell->part_idx);
        }
    }
}
