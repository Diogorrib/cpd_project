#include <math.h>
#include "simulation.h"

/**
 * @brief Get correct x index & correct x distance (wrap around)
 */
double fix_x_position_wrap_around(long *x)
{
    double dx = 0;
    if (*x < 0) {
        *x = ncside - 1;
        dx = -side;
    } else if (*x >= ncside) {
        *x = 0;
        dx = side;
    }
    return dx;
}

/**
 * @brief Get correct y distance (wrap around)
 */
double fix_y_position_wrap_around(long y)
{
    double dy = 0;
    if (rank == 0 && y == 0) {
        dy = -side;
    } else if (rank == process_count - 1 && y == block_size + 1) {
        dy = side;
    }
    return dy;
}

/**
 * @brief Get the center of mass of the adjacent cells,
 * with fixed wrap around for x and y positions
 */
void get_adj_indexes(long long cell_idx, cell_t *cells, center_of_mass_t *adj_cells)
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
        double dx = fix_x_position_wrap_around(&adj[i][0]);
        double dy = fix_y_position_wrap_around(adj[i][1]);
        long long index = adj[i][0] + adj[i][1] * ncside;

        // fix wrap around
        cell_t *cell = &cells[index];
        center_of_mass_t *adj_cell = &adj_cells[i];
        adj_cell->x = cell->x + dx;
        adj_cell->y = cell->y + dy;
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


// low number of rows per process
void compute_forces_normal(particle_t *par, cell_t *cells, MPI_Request *r)
{
    wait_for_center_of_mass(r);

    #pragma omp parallel for schedule(dynamic, 1)
    for (long long i = 0; i < n_local_cells; i++) { // Last row is ignored by n_local_cells
        long long cell_idx = i + ncside; // Skip the first row as it is computed by another process
        cell_t *cell = &cells[cell_idx];
        if (cell->n_part == 0) continue;

        center_of_mass_t adj_cells[ADJ_CELLS];
        get_adj_indexes(cell_idx, cells, adj_cells);

        for (long long j = 0; j < cell->n_part; j++) {
            particle_t *p = &par[cell->part_idx[j]];
            if (p->m != 0) {
                compute_acc_for_part(cell, p, j, par, adj_cells);
            }
        }
    }
}

// high number of rows per process
void compute_forces_maximize(particle_t *par, cell_t *cells, MPI_Request *r)
{
    long long loop_size = n_local_cells - 4*ncside; // remove 2 first and 2 last rows
    long long remove_first_row_offset = 3*ncside;

    // loop without critical rows
    #pragma omp parallel for schedule(dynamic, 1)
    for (long long i = 0; i < loop_size; i++) {
        long long cell_idx = i + remove_first_row_offset;
        cell_t *cell = &cells[cell_idx];
        if (cell->n_part == 0) continue;

        center_of_mass_t adj_cells[ADJ_CELLS];
        get_adj_indexes(cell_idx, cells, adj_cells);

        for (long long j = 0; j < cell->n_part; j++) {
            particle_t *p = &par[cell->part_idx[j]];
            if (p->m != 0) {
                compute_acc_for_part(cell, p, j, par, adj_cells);
            }
        }
    }

    // needed for critical rows
    wait_for_center_of_mass(r);

    #pragma omp parallel for schedule(dynamic, 1)
    for (long long i = 0; i < 2*ncside; i++) {
        long long cell_idx = i + n_local_cells - ncside; // last rows
        cell_t *cell = &cells[cell_idx];
        if (cell->n_part == 0) continue;

        center_of_mass_t adj_cells[ADJ_CELLS];
        get_adj_indexes(cell_idx, cells, adj_cells);

        for (long long j = 0; j < cell->n_part; j++) {
            particle_t *p = &par[cell->part_idx[j]];
            if (p->m != 0) {
                compute_acc_for_part(cell, p, j, par, adj_cells);
            }
        }
    }

    #pragma omp parallel for schedule(dynamic, 1)
    for (long long i = 0; i < 2*ncside; i++) {
        long long cell_idx = i + ncside; // first rows
        cell_t *cell = &cells[cell_idx];
        if (cell->n_part == 0) continue;

        center_of_mass_t adj_cells[ADJ_CELLS];
        get_adj_indexes(cell_idx, cells, adj_cells);

        for (long long j = 0; j < cell->n_part; j++) {
            particle_t *p = &par[cell->part_idx[j]];
            if (p->m != 0) {
                compute_acc_for_part(cell, p, j, par, adj_cells);
            }
        }
    }    
}

/**
 * @brief Compute the forces for each particle
 *
 * @param par array of local particles
 * @param cells array of local cells accounting for 2 adjacents rows
 * @param r array of MPI requests to wait to receive and to send
 */
void compute_forces(particle_t *par, cell_t *cells, MPI_Request *r)
{
    if (block_size < 5) {
        compute_forces_normal(par, cells, r);
    } else {
        compute_forces_maximize(par, cells, r);
    }
}
