#include <stdio.h>
#include "utils.h"

void* allocate_memory(size_t n_elements, size_t element_size)
{
    void *array = malloc(n_elements * element_size);
    if (!array) {
        fprintf(stderr, "Memory allocation failed\n");
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

    omp_set_lock(&cell->cell_lock);

    long long n_part = cell->n_part;
    long long *part_idx = cell->part_idx;

    long long chunk_size = get_dynamic_chunk_size(n_part);
    if (n_part == 0) {
        part_idx = (long long *)malloc(chunk_size * sizeof(long long));
    } else if (n_part % chunk_size == 0) {
        part_idx = (long long *)realloc(part_idx, (n_part + chunk_size) * sizeof(long long));
    }
    if (!part_idx) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    part_idx[n_part] = idx;
    cell->n_part++;
    cell->part_idx = part_idx;

    omp_unset_lock(&cell->cell_lock);

    par[idx].cell_idx = cell_idx;
}

void cleanup_cells(long ncside, cell_t *cells)
{
    #pragma omp for
    for (long long i = 0; i < ncside * ncside; i++) {
        cell_t *cell = &cells[i];
        if (cell->n_part > 0) {
            free(cell->part_idx);
        }
    }
}

void init_lock_cells(long ncside, cell_t *cells)
{
    #pragma omp for
    for (long long i = 0; i < ncside * ncside; i++) {
        omp_init_lock(&cells[i].cell_lock);
    }
}

void destroy_lock_cells(long ncside, cell_t *cells)
{
    #pragma omp for
    for (long long i = 0; i < ncside * ncside; i++) {
        omp_destroy_lock(&cells[i].cell_lock);
    }
}
