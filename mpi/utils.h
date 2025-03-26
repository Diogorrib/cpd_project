#ifndef UTILS_H
#define UTILS_H

#include <stdlib.h>
#include "particle.h"

void* allocate_memory(size_t n_elements, size_t element_size);
void append_particle_index(long long idx, long long cell_idx, particle_t *par, cell_t *cells);
int part_process_space(long ncside, double side, long block_low, long block_high, particle_t *p);
long long get_dynamic_chunk_size(long long n_part);
long long get_cell_idx(double inv_cell_side, long ncside, particle_t *p);
void cleanup_cells(long ncside, cell_t *cells);

#endif // UTILS_H
