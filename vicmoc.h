/*
 * VIC-MOC - vicmoc.h - structures and defines for the program
 *
 * Copyright 2004 Mark J. Stock mstock@umich.edu
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifdef __LINUX__
#include <malloc.h>
#endif
#include <math.h>
//#include <f2c.h>
#include <time.h>

#define M_PI 3.14159265358979323846

// this can and should be changed to allow running non-square sims:
// no: eventually this will support non-square *pixels*, but for now it's broken
//#define SCALE 1.7777778
#define SCALE 1.0
//#define SCALE 0.0857

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
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define TRUE 1
#define FALSE 0
#define EPSILON 1.0e-6
#define WALL 0		// these are for boundary conditions
#define OPEN 1
#define PERIODIC 2

#define VMAXAVG 5

// in libvicmoc.c
extern int add_sharp_circular_blob(int,int,int,int,float**,float,float,float,float);
extern int add_smooth_circular_blob(int,int,int,int,float**,float,float,float,float);
extern int add_smooth_circular_blob_3d(int,int,int,int,int,int,int,float***,float,float,float,float);
extern int add_smooth_spherical_blob(int,int,int,int,int,int,float***,float,float,float,float,float);
extern int add_singular_blob_3d(int,int,int,int,int,int,float***,float,float,float,float,float);
extern int add_cube_3d(int,int,int,int,int,int,float***,float,float,float,float,float);
extern int add_vortex_ring_3d(int,int,int,int,int,int,float***,float***,float***,float,float,float,float,float,float,float,float,float);
extern int write_output(char*,int,int,float**,float,float,int);
extern int write_png (char*, int, int, int, int, float**, float, float, float**, float, float, float**, float, float);
extern int read_png (char*, int, int, int, int, float, int, float**, float, float, float**, float, float, float**, float, float);
extern int write_output_3d(char*,int,int,int,float***,float,float,int,int);
extern int write_output_particles_rad(char*,int,float**,float**,float*);
extern int explicit_particle_move_3d(int,int,int,int,int,int,float****,float,float,int,float**,float**);
extern float step_forward_2d(int,int,int,int,int,int,int,int,float*,int,int,float***,float***,float***,int,float**,float,float*,int,float*,float,int,float,float,float***);
extern float step_forward_3d(int,int,int,int,int,int,int,int,float****,float****,float****,int,int*,float*,float,float);
extern int make_solenoidal_3d(int,int,int,int,int,int,float***,float***,float***);
extern float find_energy_3d(int,int,int,int,int,int,float***,float***,float***);
extern float find_vmax(float***, float***, float***,int,int,int);

extern float* allocate_1d_array_f(long int);
extern float** allocate_2d_array_f(int,int);
extern int free_2d_array_f(float**);
extern float*** allocate_3d_array_f(int,int,int);
extern int free_3d_array_f(float***);
extern float** flatten_to_2d (float***, int, int, int, int);

extern void update_mask_with_blocks (float**, float***, int, int, float);

// in gr23.c
//extern int hwscrt_(real *a, real *b, integer *m, integer *mbdcnd,
//                   real *bda, real *bdb,
//                   real *c__, real *d__, integer *n, integer *nbdcnd,
//                   real *bdc, real *bdd,
//                   real *elmbda, real *f, integer *idimf, real *pertrb,
//                   integer *ierror, real *w);

//extern int hw3crt_(real *xs, real *xf, integer *l, integer *lbdcnd, real *bdxs, real *bdxf,
//                   real *ys, real *yf, integer *m, integer *mbdcnd, real *bdys, real *bdyf,
//                   real *zs, real *zf, integer *n, integer *nbdcnd, real *bdzs, real *bdzf,
//                   real *elmbda, integer *ldimf, integer *mdimf, real *f, real *pertrb,
//                   integer *ierror, real *w);

// in mud3sp_full.f
extern
  void mud3sp_(int *iparm, float *fparm, float *work, char *cfx, char* cfy, char* cfz,
               char *bndyc, float *rhs, float *phi, int *mgopt, int *ierror);

// in mud2sp_full.f
extern
  void mud2sp_(int *iparm, float *fparm, float *work, char *cfx, char* cfy,
               char *bndyc, float *rhs, float *phi, int *mgopt, int *ierror);

// in gr2.f
extern
  void hwscrt_(float *a, float *b, int *m, int *mbdcnd,
               float *bda, float *bdb,
               float *c__, float *d__, int *n, int *nbdcnd,
               float *bdc, float *bdd,
               float *elmbda, float *f, int *idimf, float *pertrb,
               int *ierror, float *w);

//extern void MAIN__() {}

