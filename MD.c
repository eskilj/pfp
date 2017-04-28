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

        /* calculate distance from central mass */
        for (k = 0; k < Nbody; k++) {
            r[k] = 0.0;
        }

        for (i = 0; i < Nbody; i++) {
            r[i] = (pos[i][0]*pos[i][0]) + (pos[i][1]*pos[i][1]) + (pos[i][2]*pos[i][2]);
            r[i] = sqrt(r[i]);

            /* calculate central force */
            for (j = 0; j < Ndim; j++) {
                f[i][j] = - visc[i]*vel[i][j] - visc[i]*wind[j] - force(G * mass[i] * M_central, pos[i][j], r[i]);
            }
        }


        /* calculate pairwise separation of particles */
        k = 0;
        for (i = 0; i < Nbody; i++) {
            for (j = i + 1; j < Nbody; j++) {
                for (l = 0; l < Ndim; l++) {
                    delta_pos[k][l] = pos[i][l] - pos[j][l];
                }
                k = k + 1;
            }
        }

        /* calculate norm of seperation vector */
        for (k = 0; k < Npair; k++) {
            delta_r[k] = 0.0;
        }

        for (i = 0; i < Npair; i++) {
            for (j = 0; j < Ndim; j++) {
                delta_r[i] += delta_pos[i][j]*delta_pos[i][j];
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
            outer_mass = G*mass[i];
            for (j = i + 1; j < Nbody; j++) {
                inner_mass = outer_mass * mass[j];
                for (l = 0; l < Ndim; l++) {
                    /*  flip force if close in */
                    temp_force = force(inner_mass, delta_pos[k][l], delta_r[k]);
                    if (delta_r[k] >= size) {
                        f[i][l] = f[i][l] - temp_force;
                        f[j][l] = f[j][l] + temp_force;
                    } else {
                        f[i][l] = f[i][l] + temp_force;
                        f[j][l] = f[j][l] - temp_force;
                        collisions++;
                    }
                }
                k = k + 1;
            }
        }

        /* update positions
         * update velocities
         * 4473061.sdb
         */

        for (i = 0; i < Nbody; i++) {
            inner_mass = mass[i];
            for (j = 0; j < Ndim; j++) {
                pos[i][j] = pos[i][j] + dt * vel[i][j];
                vel[i][j] = vel[i][j] + dt * (f[i][j] / inner_mass);
            }
        }


    }

}




