#include <math.h>
#include "simulation.h"

void get_adj_indexes(double side, long ncside, long long cell_idx, cell_t *cells, center_of_mass_t *adj_cells)
{
    long x = cell_idx % ncside;
    long y = cell_idx / ncside;

    long adj[ADJ_CELLS][2] = {
        {x - 1, y - 1},
        {x - 1, y},
        {x - 1, y + 1},
        {x, y - 1},
        {x, y + 1},
        {x + 1, y - 1},
        {x + 1, y},
        {x + 1, y + 1}
    };

    for (int i = 0; i < ADJ_CELLS; i++) {
        double d[2] = {0, 0};
        for (int j = 0; j < 2; j++) {
            if (adj[i][j] < 0) {
                adj[i][j] = ncside - 1;
                d[j] = -side;
            } else if (adj[i][j] >= ncside) {
                adj[i][j] = 0;
                d[j] = side;
            }
        }
        long long index = adj[i][0] + adj[i][1] * ncside;
        // fix wrap around
        cell_t *cell = &cells[index];
        center_of_mass_t *adj_cell = &adj_cells[i];
        adj_cell->x = cell->x + d[0];
        adj_cell->y = cell->y + d[1];
        adj_cell->m = cell->m;
    }
}

void compute_force(double x1, double y1, double m1, double x2, double y2, double m2, double *fx, double *fy)
{
    double dx = x2 - x1;
    double dy = y2 - y1;
    double dist_sq = dx * dx + dy * dy;
    double f = G * m1 * m2 / dist_sq;

    double inv_dist = 1.0 / sqrt(dist_sq);
    double vecx = dx * inv_dist;
    double vecy = dy * inv_dist;

    *fx = f * vecx;
    *fy = f * vecy;
}

/**
 * @brief Compute the force between two particles,
 * accumulating the force in both particles (avoid duplicate calculations)
 */
void compute_force_particle_particle(particle_t *p1, particle_t *p2)
{
    double fx;
    double fy;
    compute_force(p1->x, p1->y, p1->m, p2->x, p2->y, p2->m, &fx, &fy);

    p1->fx += fx;
    p1->fy += fy;
    p2->fx -= fx;
    p2->fy -= fy;
}

/**
 * @brief Compute the force between a particle and the center of mass of a cell,
 * accumulating the force in the particle
 */
void compute_force_particle_cm(particle_t *p, center_of_mass_t *cm)
{
    double fx;
    double fy;
    compute_force(p->x, p->y, p->m, cm->x, cm->y, cm->m, &fx, &fy);

    p->fx += fx;
    p->fy += fy;
}

void compute_acc_for_part(cell_t *cell, particle_t *p1, long long j, particle_t *par, center_of_mass_t *adj_cells)
{
    // particles in the same cell
    for (long long i = j + 1; i < cell->n_part; i++) {
        particle_t *p2 = &par[cell->part_idx[i]];
        if (p2->m != 0) {
            compute_force_particle_particle(p1, p2);
        }
    }

    // center of mass of adjacent cells
    for (int i = 0; i < ADJ_CELLS; i++) {
        center_of_mass_t *adj_cell = &adj_cells[i];
        if (adj_cell->m != 0) {
            compute_force_particle_cm(p1, adj_cell);
        }
    }
}

/**
 * @brief Compute the forces for each particle and update the acceleration
 * 
 * @param side size of the side of the squared space of simulation
 * @param ncside number of cells on each side
 * @param n_part number of particles
 * @param par array of particles
 * @param cells array of cells
 */
void compute_forces(double side, long ncside, particle_t *par, cell_t *cells)
{   
    int size;
    
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    long long n_cells = ncside * ncside;
    for (long long i = 0; i < n_cells; i++) {
        cell_t *cell = &cells[i];
        if (cell->n_part == 0) continue;

        center_of_mass_t adj_cells[ADJ_CELLS];
        get_adj_indexes(side, ncside, i, cells, adj_cells);
        for (long long j = 0; j < cell->n_part; j++) {
            particle_t *p = &par[cell->part_idx[j]];
            if (p->m != 0) {
                compute_acc_for_part(cell, p, j, par, adj_cells);
            }
        }
    }
}
