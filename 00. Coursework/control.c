/*
 * $Id: control-c.c,v 1.2 2002/01/08 12:32:48 spb Exp spb $
 *
 * Control program for the MD update
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define DECL

#include "coord.h"

#define Nstep 100
#define Nsave 5
#define dt 0.02

double second(void);

int main(int argc, char *argv[]) {
    int i, j;
    FILE *in, *out;
    double tstart, tstop;
    double start, stop;
    char name[80];

    wind[Xcoord] = 0.9;
    wind[Ycoord] = 0.4;
    wind[Zcoord] = 0.0;
    /* set up multi dimensional arrays */


/* read the initial data from a file */

    collisions = 0;
    in = fopen("input.dat", "r");

    if (!in) {
        perror("input.dat");
        exit(1);
    }

    for (i = 0; i < Nbody; i++) {
        fscanf(in, "%16le%16le%16le%16le%16le%16le%16le%16le%16le\n",
               mass + i, radius + i, visc + i,
               &pos[Xcoord][i], &pos[Ycoord][i], &pos[Zcoord][i],
               &vel[Xcoord][i], &vel[Ycoord][i], &vel[Zcoord][i]);
    }
    fclose(in);

/*
 * Run Nstep timesteps and time how long it took
 */

    tstart = second();
    for (j = 1; j <= Nsave; j++) {
        start = second();
        evolve(Nstep, dt);
        stop = second();
        printf("%d timesteps took %f seconds\n", Nstep, stop - start);
        printf("collisions %d\n", collisions);

/* write final result to a file */
        sprintf(name, "output.dat%03d", j * Nstep);
        out = fopen(name, "w");

        if (!out) {
            perror(name);
            exit(1);
        }

        for (i = 0; i < Nbody; i++) {
            fprintf(out,
                    "%16.8E%16.8E%16.8E%16.8E%16.8E%16.8E%16.8E%16.8E%16.8E\n",
                    mass[i], radius[i], visc[i],
                    pos[Xcoord][i], pos[Ycoord][i], pos[Zcoord][i],
                    vel[Xcoord][i], vel[Ycoord][i], vel[Zcoord][i]);
        }
        fclose(out);
    }
    tstop = second();
    printf("%d timesteps took %f seconds\n", Nsave * Nstep, tstop - tstart);

}

double second() {
    struct timeval tp;
    struct timezone tzp;
    int i;

    i = gettimeofday(&tp, &tzp);
    return ((double) tp.tv_sec + (double) tp.tv_usec * 1.e-6);
}

