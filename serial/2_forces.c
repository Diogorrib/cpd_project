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
        adj_cells[i].x = cells[index].x + d[0];
        adj_cells[i].y = cells[index].y + d[1];
        adj_cells[i].m = cells[index].m;
    }
}

void compute_force(double x1, double y1, double m1, double x2, double y2, double m2, double *fx, double *fy)
{
    double dx = x2 - x1;
    double dy = y2 - y1;
    double dist_sq = dx * dx + dy * dy;
    double f = G * m1 * m2 / dist_sq;

    double dist = sqrt(dist_sq);
    double vecx = dx / dist;
    double vecy = dy / dist;

    *fx += f * vecx;
    *fy += f * vecy;
}

void compute_acc_for_part(double side, long ncside, long long cell_idx, long long p1_idx, particle_t *par, cell_t *cells)
{
    double fx = 0;
    double fy = 0;

    // particles in the same cell
    for (long long i = 0; i < cells[cell_idx].n_part; i++) {
        long long p2_idx = cells[cell_idx].part_idx[i];
        if (p2_idx != p1_idx && par[p2_idx].m != 0) {
            compute_force(
                par[p1_idx].x, par[p1_idx].y, par[p1_idx].m,
                par[p2_idx].x, par[p2_idx].y, par[p2_idx].m,
                &fx, &fy
            );
        }
    }

    // center of mass of adjacent cells
    center_of_mass_t adj_cells[ADJ_CELLS];
    get_adj_indexes(side, ncside, cell_idx, cells, adj_cells);
    for (int i = 0; i < ADJ_CELLS; i++) {
        if (adj_cells[i].m != 0) {
            compute_force(
                par[p1_idx].x, par[p1_idx].y, par[p1_idx].m,
                adj_cells[i].x, adj_cells[i].y, adj_cells[i].m,
                &fx, &fy
            );
        }
    }

    par[p1_idx].ax = fx / par[p1_idx].m;
    par[p1_idx].ay = fy / par[p1_idx].m;
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
    long long n_cells = ncside * ncside;
    for (long long i = 0; i < n_cells; i++) {
        for (long long j = 0; j < cells[i].n_part; j++) {
            long long p_idx = cells[i].part_idx[j];
            if (par[p_idx].m != 0) {
                compute_acc_for_part(side, ncside, i, p_idx, par, cells);
            }
        }
    }
}

/* void compute_forces(double side, long ncside, long long n_part, particle_t *par, cell_t *cells)
{
    for (long long i = 0; i < n_part; i++) {
        if (par[i].m == 0) {
            continue;
        }

        double fx = 0;
        double fy = 0;

        long long cell_idx = par[i].cell_idx;
        for (long long k = 0; k < cells[cell_idx].n_part; k++) {
            long long j = cells[cell_idx].part_idx[k];
            if (i != j && par[j].m != 0) {
                compute_force(
                    par[i].x, par[i].y, par[i].m,
                    par[j].x, par[j].y, par[j].m,
                    &fx, &fy
                );
            }
        }

        center_of_mass_t adj_cells[ADJ_CELLS];
        get_adj_indexes(side, ncside, cell_idx, cells, adj_cells);
        for (int k = 0; k < ADJ_CELLS; k++) {
            if (adj_cells[k].m != 0) {
                compute_force(
                    par[i].x, par[i].y, par[i].m,
                    adj_cells[k].x, adj_cells[k].y, adj_cells[k].m,
                    &fx, &fy
                );
            }
        }

        par[i].ax = fx / par[i].m;
        par[i].ay = fy / par[i].m;
    }
} */
