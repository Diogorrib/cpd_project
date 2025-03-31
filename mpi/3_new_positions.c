#include "simulation.h"
#include <stdio.h>

/**
 * @brief Distributes particles into cells based on their positions,
 * excluding particles that have collided (mass = 0)
 *
 * @param par array of particles
 * @param cells array of cells
 */
void particle_distribution(particle_t *par, cell_t *cells)
{
    // init / reset cells
    for (long long i = 0; i < n_local_cells; i++) { // Last row is ignored by n_local_cells
        long long cell_idx = i + ncside; // Skip the first row as it is computed by another process
        cells[cell_idx].n_part = 0;
    }

    for (long long i = 0; i < n_part; i++) {
        particle_t *p = &par[i];
        if (p->m == 0) continue;

        long long cell_idx = get_local_cell_idx(p);

        if(cell_in_process_space(cell_idx)) {
            append_particle_index(i, cell_idx, par, cells);
        } else {
        }
    }
    //TODO: Send particles that moved to another process(ignore in the first execution)
}

void check_outside_space(particle_t *p)
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
 * @param par array of particles
 * @param cells array of cells
 */
void compute_new_positions(particle_t *par, cell_t *cells)
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
        check_outside_space(p);

        // reset force for next iteration
        p->fx = 0;
        p->fy = 0;
    }
    cleanup_cells(cells);
    particle_distribution(par, cells);
}
