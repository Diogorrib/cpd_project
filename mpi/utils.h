#ifndef UTILS_H
#define UTILS_H

#include <stdlib.h>
#include "particle.h"

#define BLOCK_LOW(id, p, ncside) ((id) * (ncside) / (p))
#define BLOCK_HIGH(id, p, ncside) (BLOCK_LOW((id) + 1, p, ncside) - 1)
#define BLOCK_SIZE(id, p, ncside) (BLOCK_HIGH(id, p, ncside) - BLOCK_LOW(id, p, ncside) + 1)
#define BLOCK_OWNER(index, p, ncside) (((p) * ((index) + 1) - 1) / (ncside))

void* allocate_memory(size_t n_elements, size_t element_size);
void append_particle_to_array(long long idx, particle_t *p, particle_t **p_array);
void append_particle_to_cell(long long idx, long long cell_idx, particle_t *par, cell_t *cells);
int cell_in_process_space(long long cell_idx);
long long get_local_cell_idx(particle_t *p);
long long get_dynamic_chunk_size(long long n_part);
void cleanup_cells(cell_t *cells);

#endif // UTILS_H
