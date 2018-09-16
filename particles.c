/*
 * VIC-MOC - particles.c - create, move, draw particles
 *
 * Copyright 2018 Mark J. Stock mstock@umich.edu
 */

#include <stdlib.h>
#include <stdio.h>
#include "particles.h"
#include "utility.h"

// external functions
extern int interpolate_vel_using_M4p_2d (int,int,int,int,float**,float**,float**,float,float,float*,float*);

// internal functions
int resize_particles (struct Particles*, const int);

/*
 * initialize the struct
 */
int init_particles (struct Particles *p, const int nspace) {
   p->n = 0;
   p->nmax = nspace;
   p->x = allocate_1d_array_f(2*nspace);
   p->u = allocate_1d_array_f(2*nspace);
   p->c = allocate_1d_array_f(3*nspace);
   p->m = allocate_1d_array_f(nspace);
   p->b = allocate_1d_array_f(nspace);
   //(void) resize_particles(p, nspace);
   return p->n;
}

/*
 * resize the struct
 */
int resize_particles (struct Particles *p, const int nspace) {

   // need to save the previous values!!!
   float* temp = p->x;
   p->x = allocate_1d_array_f(2*nspace);
   for (int i=0; i<2*p->nmax; ++i) p->x[i] = temp[i];
   free_1d_array_f(temp);

   temp = p->u;
   p->u = allocate_1d_array_f(2*nspace);
   for (int i=0; i<2*p->nmax; ++i) p->u[i] = temp[i];
   free_1d_array_f(temp);

   temp = p->c;
   p->c = allocate_1d_array_f(3*nspace);
   for (int i=0; i<3*p->nmax; ++i) p->c[i] = temp[i];
   free_1d_array_f(temp);

   temp = p->m;
   p->m = allocate_1d_array_f(nspace);
   for (int i=0; i<p->nmax; ++i) p->m[i] = temp[i];
   free_1d_array_f(temp);

   temp = p->b;
   p->b = allocate_1d_array_f(nspace);
   for (int i=0; i<p->nmax; ++i) p->b[i] = temp[i];
   free_1d_array_f(temp);

   p->nmax = nspace;
   return 0;
}

/*
 * generate a block of particles
 */
int add_block_of_particles (struct Particles *p, int nnew,
                            float xs, float xf, float ys, float yf,
                            float rr, float gg, float bb,
                            float m, float b) {

   // make sure we have room
   if (p->n+nnew > p->nmax) resize_particles(p, p->n+nnew);

   // fill in the particles!
   for (int i=p->n; i<p->n+nnew; ++i) {
      // position and vel
      p->x[2*i+0] = xs + (xf-xs)*(rand()/(float)RAND_MAX);
      p->x[2*i+1] = ys + (yf-ys)*(rand()/(float)RAND_MAX);
      p->u[2*i+0] = 0.0;
      p->u[2*i+1] = 0.0;
      // color
      p->c[3*i+0] = rr;
      p->c[3*i+1] = gg;
      p->c[3*i+2] = bb;
      // dynamics parameters
      p->m[i] = m;
      p->b[i] = b;
   }

   // update the active count
   p->n += nnew;

   return p->n;
}

/*
 * move the particles along the local flow, use 2nd order backward method
 */
int move_particles (struct Particles *p,
                    int nx, int ny, int xbdry, int ybdry,
                    float **mask, float **u, float **v, float **temp,
                    float *grav, const float dt) {

   float spx, spy, newx, newy;
   float u0, v0, u1, v1;

   for (int i=0; i<p->n; ++i) {
      spx = p->x[2*i+0];
      spy = p->x[2*i+1];
      interpolate_vel_using_M4p_2d(nx,ny,xbdry,ybdry,mask,u,v,spx,spy,&u0,&v0);
      // must somehow apply effects of momentum force (m) and buoyancy force (b)
      newx = spx - dt*u0;
      newy = spy - dt*v0;
      interpolate_vel_using_M4p_2d(nx,ny,xbdry,ybdry,mask,u,v,newx,newy,&u1,&v1);
      // compose new velocity as weighted average of original and fluid velocities (do we need dt in this?)
      p->u[2*i+0] = (p->m[i]*p->u[2*i+0] + 0.5*(u0+u1)) / (1.0+p->m[i]);
      p->u[2*i+1] = (p->m[i]*p->u[2*i+0] + 0.5*(v0+v1)) / (1.0+p->m[i]);
   }

   // apply the new velocity
   for (int i=0; i<2*p->n; ++i) {
      p->x[i] += dt*p->u[i];
   }

   //fprintf(stdout,"  particle 1 is at %g %g with velocity %g %g\n",p->x[0],p->x[1],p->u[0],p->u[1]);

   return p->n;
}

/*
 * splat particles onto a color image
 */
int draw_particles (struct Particles *p, float yf, int nx, int ny,
                    float **inred, float **ingrn, float **inblu,
                    float **outred, float **outgrn, float **outblu) {

   // copy the in colors to the out array
   for (int ix=0; ix<nx; ix++) {
      for (int iy=0; iy<ny; iy++) {
         outred[ix][iy] = inred[ix][iy];
         outgrn[ix][iy] = ingrn[ix][iy];
         outblu[ix][iy] = inblu[ix][iy];
      }
   }

   // splat the particles over this
   //fprintf(stdout,"  splatting %d particles\n",p->n);
   for (int i=0; i<p->n; ++i) {
      const float px = p->x[2*i+0] * (float)(nx-1);
      const float py = p->x[2*i+1] * (float)(ny-1) / yf;
      const int ix = (int)px;
      const int iy = (int)py;
      const float fx = px - (float)ix;
      const float fy = py - (float)iy;
      // drop the color on weighted by the ballistic coefficient (so heavy particles draw more heavily)
      const float wgt = p->m[i];
      //if (i==0) fprintf(stdout,"  particle 1 is at %g %g with index %d %d and weight %g\n",px,py,ix,iy,wgt);
      float fac = (1.0-fx)*(1.0-fy);
      outred[ix][iy] = (outred[ix][iy] + fac*wgt*p->c[3*i+0]) / (1.0+fac*wgt);
      outgrn[ix][iy] = (outgrn[ix][iy] + fac*wgt*p->c[3*i+1]) / (1.0+fac*wgt);
      outblu[ix][iy] = (outblu[ix][iy] + fac*wgt*p->c[3*i+2]) / (1.0+fac*wgt);
      fac = fx*(1.0-fy);
      outred[ix+1][iy] = (outred[ix+1][iy] + fac*wgt*p->c[3*i+0]) / (1.0+fac*wgt);
      outgrn[ix+1][iy] = (outgrn[ix+1][iy] + fac*wgt*p->c[3*i+1]) / (1.0+fac*wgt);
      outblu[ix+1][iy] = (outblu[ix+1][iy] + fac*wgt*p->c[3*i+2]) / (1.0+fac*wgt);
      fac = fy*(1.0-fx);
      outred[ix][iy+1] = (outred[ix][iy+1] + fac*wgt*p->c[3*i+0]) / (1.0+fac*wgt);
      outgrn[ix][iy+1] = (outgrn[ix][iy+1] + fac*wgt*p->c[3*i+1]) / (1.0+fac*wgt);
      outblu[ix][iy+1] = (outblu[ix][iy+1] + fac*wgt*p->c[3*i+2]) / (1.0+fac*wgt);
      fac = fy*fx;
      outred[ix+1][iy+1] = (outred[ix+1][iy+1] + fac*wgt*p->c[3*i+0]) / (1.0+fac*wgt);
      outgrn[ix+1][iy+1] = (outgrn[ix+1][iy+1] + fac*wgt*p->c[3*i+1]) / (1.0+fac*wgt);
      outblu[ix+1][iy+1] = (outblu[ix+1][iy+1] + fac*wgt*p->c[3*i+2]) / (1.0+fac*wgt);
      //if (iy > -1) {
      //}
   }
}

