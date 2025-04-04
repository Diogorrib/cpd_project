#ifndef GLOBALS_H
#define GLOBALS_H

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

// index of the particle zero
extern long long particle_0_idx;

extern int rank, process_count;
extern long block_low, block_size;
extern int prev_rank, next_rank;

#endif // GLOBALS_H
