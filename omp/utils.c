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

void append_particle_index(long long idx, long long cell_idx, particle_t *par, cell_t *cells)
{
    cell_t *cell = &cells[cell_idx];

    omp_set_lock(&cell->cell_lock);

    long long n_part = cell->n_part;
    long long *part_idx = cell->part_idx;

    if (n_part == 0) {
        part_idx = (long long *)malloc(sizeof(long long));
    } else {
        part_idx = (long long *)realloc(part_idx, (n_part + 1) * sizeof(long long));
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
    #pragma omp parallel for
    for (long long i = 0; i < ncside * ncside; i++) {
        cell_t *cell = &cells[i];
        if (cell->n_part > 0) {
            free(cell->part_idx);
        }
    }
}
