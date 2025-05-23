#define _USE_MATH_DEFINES
#include <math.h>
#include "utils.h"
#include "constant.h"
#include "particle.h"

unsigned int seed;
void init_r4uni(int input_seed)
{
    seed = input_seed + 987654321;
}
double rnd_uniform01()
{
    int seed_in = seed;
    seed ^= (seed << 13);
    seed ^= (seed >> 17);
    seed ^= (seed << 5);
    return 0.5 + 0.2328306e-09 * (seed_in + (int) seed);
}
double rnd_normal01()
{
    double u1, u2, z, result;
    do {
        u1 = rnd_uniform01();
        u2 = rnd_uniform01();
        z = sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
        result = 0.5 + 0.15 * z;      // Shift mean to 0.5 and scale
    } while (result < 0 || result >= 1);
    return result;
}

long long init_particles(long userseed, double side, long ncside, long long n_part, long long *part_0_idx, particle_t **par)
{
    double (*rnd01)() = rnd_uniform01;
    long long i, j = 0;

    if(userseed < 0) {
        rnd01 = rnd_normal01;
        userseed = -userseed;
    }
    
    init_r4uni(userseed);

    for(i = 0; i < n_part; i++) {
        particle_t p;
        p.x = rnd01() * side;
        p.y = rnd01() * side;
        p.vx = (rnd01() - 0.5) * side / ncside / 5.0;
        p.vy = (rnd01() - 0.5) * side / ncside / 5.0;

        p.m = rnd01() * 0.01 * (ncside * ncside) / n_part / G * EPSILON2;
        long long cell_idx = get_local_cell_idx(&p);
        if(cell_in_process_space(cell_idx)) {
            p.is_particle_0 = (i == 0) ? 1 : 0;

            append_particle_to_array(j, &p, par, 0);

            j++;
        }
    }

    // save index of the first particle
    if(j > 0 && (*par)[0].is_particle_0) {
        *part_0_idx = 0;
    }
    return j;
}
