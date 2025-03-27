#include "simulation.h"
#include <stdio.h>

/**
 * @brief Distributes particles into cells based on their positions,
 * excluding particles that have collided (mass = 0)
 * 
 * @param side size of the side of the squared space of simulation
 * @param ncside number of cells on each side
 * @param block_size number of cells in the block
 * @param n_part number of particles
 * @param par array of particles
 * @param cells array of cells
 */
void particle_distribution(double side, long ncside, long long block_size, long block_low, long long n_part, particle_t *par, cell_t *cells)
{
    double inv_cell_side = ncside / side;
    long long n_cells = (ncside+2) * block_size;

    // init / reset cells
    for (long long i = 0; i < n_cells; i++) {
        cells[i].n_part = 0;
    }

    for (long long i = 0; i < n_part; i++) {
        particle_t *p = &par[i];
        if (p->m == 0) continue;

        long long cell_idx = get_cell_idx(inv_cell_side, ncside, p);
        
        if(cell_process_space(ncside, block_low, block_size, cell_idx)){
            cell_idx += -block_low*ncside + ncside;

            append_particle_index(i, cell_idx, par, cells);
        } else {
        }
    }
    //TODO: Send particles that moved to another process(ignore in the first execution)
}

void check_outside_space(double side, particle_t *p)
{
    if (p->x < 0) {
        p->x += side;
    } else if (p->x >= side) {
        p->x -= side;
    }

    if (p->y < 0) {
        p->y += side;
    } else if (p->y >= side) {
        p->y -= side;
    }
}

/**
 * @brief  Compute the new velocity and position for each particle,
 * check if new positions are outside the space and update the cells
 * 
 * @param side size of the side of the squared space of simulation
 * @param ncside number of cells on each side
 * @param block_size number of cells in the block
 * @param n_part number of particles
 * @param par array of particles
 * @param cells array of cells
 */
void compute_new_positions(double side, long ncside, long long block_size, long block_low, long long n_part, particle_t *par, cell_t *cells)
{
    for (long long i = 0; i < n_part; i++) {
        particle_t *p = &par[i];
        if (p->m == 0) continue;

        double inv_mass = 1.0 / p->m;
        double ax = p->fx * inv_mass;
        double ay = p->fy * inv_mass;

        p->x += p->vx * DELTAT + 0.5 * ax * DELTAT * DELTAT;
        p->y += p->vy * DELTAT + 0.5 * ay * DELTAT * DELTAT;
        p->vx += ax * DELTAT;
        p->vy += ay * DELTAT;
        check_outside_space(side, p);

        // reset force for next iteration
        p->fx = 0;
        p->fy = 0;
    }
    cleanup_cells(ncside, block_size, cells);
    particle_distribution(side, ncside, block_size, block_low, n_part, par, cells);
}
