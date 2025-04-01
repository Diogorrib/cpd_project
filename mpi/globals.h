#ifndef GLOBALS_H
#define GLOBALS_H

#include "particle.h"

// size of the side of the squared space of simulation
extern double side;

// size of the grid (number of cells on each side)
extern long ncside;

// inversed size of the cell
extern double inv_cell_side;

// local number of particles
extern long long n_part;

// local number of cells
extern long long n_local_cells;

// reference to the particle zero
extern particle_t *particle_0;

extern int rank, process_count;
extern long block_low, block_size;//, block_high, block_size;

#endif // GLOBALS_H
