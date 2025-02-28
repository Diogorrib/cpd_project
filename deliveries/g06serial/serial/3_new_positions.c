#include "simulation.h"

/**
 * @brief Distributes particles into cells based on their positions,
 * excluding particles that have collided (mass = 0)
 * 
 * @param side size of the side of the squared space of simulation
 * @param ncside number of cells on each side
 * @param n_part number of particles
 * @param par array of particles
 * @param cells array of cells
 */
void particle_distribution(double side, long ncside, long long n_part, particle_t *par, cell_t *cells)
{
    double inv_cell_side = ncside / side;
    long long n_cells = ncside * ncside;

    // init / reset cells
    for (long long i = 0; i < n_cells; i++) {
        cells[i].n_part = 0;
    }

    for (long long i = 0; i < n_part; i++) {
        particle_t *p = &par[i];
        if (p->m == 0) continue;

        long cell_x_idx = (long)(p->x * inv_cell_side);
        long cell_y_idx = (long)(p->y * inv_cell_side);
        long long cell_idx = cell_x_idx + cell_y_idx * ncside;

        append_particle_index(i, cell_idx, par, cells);
    }
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
 * @param n_part number of particles
 * @param par array of particles
 * @param cells array of cells
 */
void compute_new_positions(double side, long ncside, long long n_part, particle_t *par, cell_t *cells)
{
    for (long long i = 0; i < n_part; i++) {
        particle_t *p = &par[i];
        if (p->m == 0) continue;

        p->x += p->vx * DELTAT + 0.5 * p->ax * DELTAT * DELTAT;
        p->y += p->vy * DELTAT + 0.5 * p->ay * DELTAT * DELTAT;
        p->vx += p->ax * DELTAT;
        p->vy += p->ay * DELTAT;
        check_outside_space(side, p);
    }
    cleanup_cells(ncside, cells);
    particle_distribution(side, ncside, n_part, par, cells);
}
