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

double force(double W, double delta, double r);

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

        for (i = 0; i < Nbody; i++) {
            /* calculate distance from central mass */
            r[i] = (pos[i][0]*pos[i][0]) + (pos[i][1]*pos[i][1]) + (pos[i][2]*pos[i][2]);
            r[i] = sqrt(r[i]);

            /* calculate central force */
            for (j = 0; j < Ndim; j++) {
                f[i][j] = -visc[i]*(vel[i][j] + wind[j]) - force(G * mass[i] * M_central, pos[i][j], r[i]);
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
        for (i = 0; i < Npair; i++) {
            delta_r[i] = (delta_pos[i][0]*delta_pos[i][0]) + (delta_pos[i][1]*delta_pos[i][1]) + (delta_pos[i][2]*delta_pos[i][2]);
            delta_r[i] = sqrt(delta_r[i]);
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




