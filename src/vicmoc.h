/*
 * VIC-MOC - vicmoc.h - structures and defines for both codes
 *
 * Copyright 2004-7,20,23 Mark J. Stock <mstock@umich.edu>
 */

#pragma once

#include <stdbool.h>

#define M_PI 3.14159265358979323846
#define MAXCHARS 255

// these can be changed
// for ordering of the arrays in "a", they will change depending on what is being tracked
#define XV 0		// x-velocity
#define YV 1		// y-velocity
#define ZV 2		// z-velocity
#define MAX_SCALARS 12	// must be greater than 3 for vic3d->moc to work
#define WX 0		// vorticity, x-dir, 3D field
#define WY 1		// vorticity, y-dir, 3D field
#define WZ 2		// vorticity, z-dir, 3D field
#define W2 0		// vorticity, 2D field
#define SF 3		// scalar fraction (density)
#define DF 4		// dilatation field
#define TF 5		// temperature field
#define RR 6		// red
#define GG 7		// green
#define BB 8		// blue
#define MD 9		// viscosity (momentum diffusivity) field (for variable viscosity)
#define TD 10		// scalar diffusivity field (for variable viscosity)
#define CD 11		// color diffusivity field (for variable viscosity)

// these don't need to be modified
#define FLOAT float	// use 4-byte floats or 8-byte doubles; actually, this must
 			// be set to "float" because the gr2.c routines require it

// see https://stackoverflow.com/questions/3437404/min-and-max-in-c for a discussion of this
#define MAX(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b;       \
})

#define MIN(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b;       \
})

#define EPSILON 1.0e-6

// these are for boundary conditions - move these to an enum
#define WALL 0
#define OPEN 1
#define PERIODIC 2

#define VMAXAVG 5

// available interpolation methods
typedef enum interpMeth {
   cic,		// 1st order cloud-in-cell (bi- or tri-linear interpolation), support width=2
   tsc,		// 2nd order triangle-shaped cloud, support width=3
   m4p,		// 4th order M4', support width=4
   cic2,
   cic3
} INTERP;


// from mud2sp_extern.c
extern void funcbndyc (int*, float*, float*, float*);
extern void funccfx (float*, float*, float*, float*);
extern void funccfy (float*, float*, float*, float*);

// from, well, color stuff
extern void get_random_color (float***, int, int, float*);

