#ifndef UTILS_H
#define UTILS_H

#include <stdlib.h>
#include "particle.h"

void* allocate_memory(size_t n_elements, size_t element_size);
void append_particle_index(long long idx, long long cell_idx, particle_t *par, cell_t *cells);
void cleanup_cells(long ncside, cell_t *cells);
void init_lock_cells(long ncside, cell_t *cells);
void destroy_lock_cells(long ncside, cell_t *cells);

#endif // UTILS_H
