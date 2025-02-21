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
    double cell_side = side / ncside;
    long long n_cells = ncside * ncside;

    // init / reset cells
    for (long long i = 0; i < n_cells; i++) {
        cells[i].n_part = 0;
    }

    for (long long i = 0; i < n_part; i++) {
        if (par[i].m == 0) {
            continue;
        }

        long cell_x_idx = (long)(par[i].x / cell_side);
        long cell_y_idx = (long)(par[i].y / cell_side);
        long long cell_idx = cell_x_idx + cell_y_idx * ncside;

        append_particle_index(i, cell_idx, par, cells);
    }
}

void check_outside_space(double side, long long p_idx, particle_t *par)
{
    if (par[p_idx].x < 0) {
        par[p_idx].x += side;
    } else if (par[p_idx].x >= side) {
        par[p_idx].x -= side;
    }

    if (par[p_idx].y < 0) {
        par[p_idx].y += side;
    } else if (par[p_idx].y >= side) {
        par[p_idx].y -= side;
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
        if (par[i].m == 0) {
            continue;
        }
        par[i].x += par[i].vx * DELTAT + 0.5 * par[i].ax * DELTAT * DELTAT;
        par[i].y += par[i].vy * DELTAT + 0.5 * par[i].ay * DELTAT * DELTAT;
        par[i].vx += par[i].ax * DELTAT;
        par[i].vy += par[i].ay * DELTAT;
        check_outside_space(side, i, par);
    }
    cleanup_cells(ncside, cells);
    particle_distribution(side, ncside, n_part, par, cells);
}
