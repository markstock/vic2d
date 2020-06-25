/*
 * VIC-MOC - particles.h - create, move, draw particles
 *
 * Copyright 2018-20 Mark J. Stock mstock@umich.edu
 */

#pragma once

struct Particles {
   // how many active?
   int n;
   // how much room in arrays
   int nmax;

   // position (2*n)
   float* x;
   // last position (2*n)
   float* oldx;
   // velocity (2*n)
   float* u;
   // color (3*n)
   float* c;
   // ballistic coefficient (like mass) - resistance to momentum (n)
   float* m;
   // buoyancy - reaction to density/temperature (n)
   float* b;

   // other ideas:
   //   unstick velocity (imposed flow must be moving this fast to get particle to move)
   //   draw radius (can allow very subtle particles as well as large ones)
   //   feedback (magnitude of effect of the particle on the flow)
};

int init_particles (struct Particles*, const int);
int add_one_particle (struct Particles*, float, float, float, float, float, float, float, float, float);
int add_block_of_particles (struct Particles*, int, float, float, float, float, float, float, float, float, float);
int move_particles (struct Particles*, int, int, int, int, float, float**, float**, float**, float**, float*, const float);
int draw_particles (struct Particles*, float, int, int, float, float**, float**, float**, int, int, float, float**, float**, float**, float, float, float);
void mult_part_by_mask (float, int, int, float**, int, int, float**, float**, float**);

