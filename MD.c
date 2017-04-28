/*
 *  Simple molecular dynamics code.
 *  $Id: MD-c.c,v 1.2 2002/01/31 16:43:14 spb Exp spb $
 *
 * This program implements:
 *     long range inverse square forces between particles. F = G * m1*m2 / r**2
 *     viscosity term     F = -u V
 * If 2 particles approach closer than Size we flip the direction of the
 * interaction force to approximate a collision.
 *
 * Coordinates are relative to a large central mass and the entire system is moving relative to the
 * viscous media.
 * If 2 particles approach closer than Size we flip the direction of the
 * interaction force to approximate a collision.
 *
 * This program was developed as part of a code optimisation course
 * and is therefore deliberately inefficient.
 *
 */
#include <stdio.h>
#include <math.h>
#include "coord.h"

void visc_force(int N, double *f, double *visc, double *vel);

void add_norm(int N, double *r, double *delta);

double force(double W, double delta, double r);

void wind_force(int N, double *f, double *visc, double vel);

#define size 1.0

void evolve(int count, double dt) {
    int step;
    int i, j, k, l;
    double inner_mass, outer_mass, temp_force;
    /*
     * Loop over timesteps.
     */
    for (step = 1; step <= count; step++) {
        printf("timestep %d\n", step);
        printf("collisions %d\n", collisions);

        /* set the viscosity term in the force calculation
         * add the wind term in the force calculation
         */

        for (j = 0; j < Ndim; j++) {
            for (i = 0; i < Nbody; i++) {
                f[j][i] = -visc[i] * vel[j][i] - visc[i] * wind[j];
            }
        }


        /* calculate distance from central mass */
        for (k = 0; k < Nbody; k++) {
            r[k] = 0.0;
        }

        for (j = 0; j < Ndim; j++) {
            for (i = 0; i < Nbody; i++) {
                r[i] += pos[j][i] * pos[j][i];
            }
        }

        for (k = 0; k < Nbody; k++) {
            r[k] = sqrt(r[k]);
        }
        /* calculate central force */
        for (j = 0; j < Ndim; j++) {
            for (i = 0; i < Nbody; i++) {
                f[j][i] = f[j][i] -
                          force(G * mass[i] * M_central, pos[j][i], r[i]);
            }
        }
        /* calculate pairwise separation of particles */
        k = 0;
        for (i = 0; i < Nbody; i++) {
            for (j = i + 1; j < Nbody; j++) {
                for (l = 0; l < Ndim; l++) {
                    delta_pos[l][k] = pos[l][i] - pos[l][j];
                }
                k = k + 1;
            }
        }

        /* calculate norm of seperation vector */
        for (k = 0; k < Npair; k++) {
            delta_r[k] = 0.0;
        }

        for (j = 0; j < Ndim; j++) {
            for (i = 0; i < Npair; i++) {
                delta_r[i] += delta_pos[j][i] * delta_pos[j][i];
            }
        }

        for (k = 0; k < Npair; k++) {
            delta_r[k] = sqrt(delta_r[k]);
        }

        /*
         * add pairwise forces.
         */
        k = 0;
        for (i = 0; i < Nbody; i++) {
            outer_mass = G * mass[i];
            for (j = i + 1; j < Nbody; j++) {
                inner_mass = outer_mass * mass[j];
                for (l = 0; l < Ndim; l++) {
                    /*  flip force if close in */
                    temp_force = force(inner_mass, delta_pos[l][k],
                                       delta_r[k]);
                    if (delta_r[k] >= size) {
                        f[l][i] = f[l][i] - temp_force;
                        f[l][j] = f[l][j] + temp_force;
                    } else {
                        f[l][i] = f[l][i] + temp_force;
                        f[l][j] = f[l][j] - temp_force;
                        collisions++;
                    }
                }
                k = k + 1;
            }
        }

        /* update positions */
        for (j = 0; j < Ndim; j++) {
            for (i = 0; i < Nbody; i++) {
                pos[j][i] = pos[j][i] + dt * vel[j][i];
            }
        }

        /* update velocities */
        for (j = 0; j < Ndim; j++) {
            for (i = 0; i < Nbody; i++) {
                vel[j][i] = vel[j][i] + dt * (f[j][i] / mass[i]);
            }
        }


    }

}




