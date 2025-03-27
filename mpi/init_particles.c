#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
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

long long init_particles(long userseed, double side, long ncside, long long n_part, long block_low, long block_high, int *first, particle_t **par)
{
    double (*rnd01)() = rnd_uniform01;
    long long chunk_size, i, j = 0;

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
        if(part_process_space(ncside, side, block_low, block_high, &p)) {
            if(i == 0) {
                (*first) = 1;    
            }

            chunk_size = get_dynamic_chunk_size(j);
            
            if (j == 0) {
                *par = (particle_t *)malloc(chunk_size * sizeof(particle_t));
            } else if (j % chunk_size == 0) {
                *par = (particle_t *)realloc(*par, (j + chunk_size) * sizeof(particle_t));
            }
            if (!*par) {
                fprintf(stderr, "Memory allocation failed (1)\n");
                exit(EXIT_FAILURE);
            }
            
            (*par)[j] = p;
            j++;
        }
    }
    return j;
}
