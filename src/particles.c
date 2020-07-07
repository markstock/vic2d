/*
 * VIC-MOC - particles.c - create, move, draw particles
 *
 * Copyright 2018-20 Mark J. Stock <mstock@umich.edu>
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <signal.h>
#include "vicmoc.h"
#include "particles.h"
#include "utility.h"

// external functions
//extern float interpolate_using_M4p_2d(int,int,int,int,float**,float,float);
extern float interpolate_array_using_M4p_2d (int,int,int,int,float**,float***,float,float,int,float*);
extern float interpolate_array_using_TSC_2d (int,int,int,int,float***,float,float,int,float*);
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
   p->oldx = allocate_1d_array_f(2*nspace);
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

   temp = p->oldx;
   p->oldx = allocate_1d_array_f(2*nspace);
   for (int i=0; i<2*p->nmax; ++i) p->oldx[i] = temp[i];
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
 * add one particle
 */
int add_one_particle (struct Particles *p,
                      float x, float y, float u, float v,
                      float rr, float gg, float bb,
                      float m, float b) {

   // make sure we have room
   if (p->n+1 > p->nmax) resize_particles(p, (int)(1.2*(p->n+100)));

   int i = p->n;

   // position and vel
   p->x[2*i+0] = x;
   p->x[2*i+1] = y;
   p->oldx[2*i+0] = x;
   p->oldx[2*i+1] = y;
   p->u[2*i+0] = u;
   p->u[2*i+1] = v;
   // color
   p->c[3*i+0] = rr;
   p->c[3*i+1] = gg;
   p->c[3*i+2] = bb;
   // dynamics parameters
   p->m[i] = m;
   p->b[i] = b;

   // update the active count
   p->n++;

   return p->n;
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
      p->oldx[2*i+0] = p->x[2*i+0];
      p->oldx[2*i+1] = p->x[2*i+1];
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
                    int nx, int ny, int xbdry, int ybdry, float yf,
                    float **mask, float **u, float **v, float **temp,
                    float *grav, const float dt) {

   float spx, spy, newx, newy;
   float u0, v0, u1, v1, u2, v2;

   #pragma omp parallel for private(spx,spy,newx,newy,u0,v0,u1,v1,u2,v2)
   for (int i=0; i<p->n; ++i) {
      // find velocity right here
      spx = p->x[2*i+0];
      spy = p->x[2*i+1];
      interpolate_vel_using_M4p_2d(nx,ny,xbdry,ybdry,mask,u,v,spx,spy,&u0,&v0);
      // now find velocity by looking back a half time step
      newx = spx - 0.5*dt*u0;
      newy = spy - 0.5*dt*v0;
      interpolate_vel_using_M4p_2d(nx,ny,xbdry,ybdry,mask,u,v,newx,newy,&u1,&v1);
      // now find velocity by looking forward a half time step
      newx = spx + 0.5*dt*u0;
      newy = spy + 0.5*dt*v0;
      interpolate_vel_using_M4p_2d(nx,ny,xbdry,ybdry,mask,u,v,newx,newy,&u2,&v2);
      // compose new velocity as weighted average of forward and backward velocities
      // and using ballistic coefficient (higher means more resistant to change)
      const float fac = p->b[i]*p->m[i];
      p->u[2*i+0] = (fac*p->u[2*i+0] + 0.5*(u2+u1)) / (1.0+fac);
      p->u[2*i+1] = (fac*p->u[2*i+1] + 0.5*(v2+v1)) / (1.0+fac);
   }

   // save the old positions
   for (int i=0; i<2*p->n; ++i) {
      p->oldx[i] = p->x[i];
   }

   // apply the new velocity
   for (int i=0; i<2*p->n; ++i) {
      p->x[i] += dt*p->u[i];
   }

   // and make sure particles stay in bounds
   if (xbdry == PERIODIC) {
      #pragma omp parallel for
      for (int i=0; i<p->n; ++i) {
         if (p->x[2*i+0] < 0.0) { p->x[2*i+0] += 1.0; p->oldx[2*i+0] += 1.0; }
         if (p->x[2*i+1] < 0.0) { p->x[2*i+1] += yf;  p->oldx[2*i+1] += yf; }
         if (p->x[2*i+0] > 1.0) { p->x[2*i+0] -= 1.0; p->oldx[2*i+0] -= 1.0; }
         if (p->x[2*i+1] > yf)  { p->x[2*i+1] -= yf;  p->oldx[2*i+1] -= yf; }
      }
   } else {
      // wall or open
      #pragma omp parallel for
      for (int i=0; i<p->n; ++i) {
         if (p->x[2*i+0] < 0.0) p->x[2*i+0] = 1.e-4;
         if (p->x[2*i+1] < 0.0) p->x[2*i+1] = 1.e-4;
         if (p->x[2*i+0] > 1.0) p->x[2*i+0] = 1.0-1.e-4;
         if (p->x[2*i+1] > yf) p->x[2*i+1] = yf-1.e-4;
      }
   }

   //fprintf(stdout,"  particle 1 is at %g %g with velocity %g %g\n",p->x[0],p->x[1],p->u[0],p->u[1]);

   return p->n;
}

/*
 * splat particles onto a color image
 *
 * input color image is nx by ny
 * particle color image is pnx by pny
 */
void draw_particles (struct Particles *p, float yf,
                     int nx, int ny,
                     float infac,  float **inred,  float **ingrn,  float **inblu,
                     int pnx, int pny,
                     float outfac, float **outred, float **outgrn, float **outblu,
                     float draw_fac, float mass_pow, float vel_pow) {

   // scale the particle (output) colors first
   for (int ix=0; ix<pnx; ix++) {
      for (int iy=0; iy<pny; iy++) {
         outred[ix][iy] *= outfac;
         outgrn[ix][iy] *= outfac;
         outblu[ix][iy] *= outfac;
      }
   }

   // then copy over colors with the input array
   if (infac != 0.0) {
      float** tempin[3];
      tempin[0] = inred;
      tempin[1] = ingrn;
      tempin[2] = inblu;
      float outvals[3];
      for (int ix=0; ix<pnx; ix++) {
         for (int iy=0; iy<pny; iy++) {
            const float px = ix / (float)(pnx-1);
            const float py = yf * iy / (float)(pny-1);
            // pick color from the input map at nx by ny (use the interp array call?)
            //const float vred = interpolate_using_M4p_2d(nx,ny,PERIODIC,PERIODIC,inred,px,py);
            interpolate_array_using_M4p_2d(nx,ny,PERIODIC,PERIODIC,NULL,tempin,px,py,3,outvals);
            // apply it to the output map at pnx by pny
            outred[ix][iy] += infac*outvals[0];
            outgrn[ix][iy] += infac*outvals[1];
            outblu[ix][iy] += infac*outvals[2];
         }
      }
   }

   // splat the particles over this
   //fprintf(stdout,"  splatting %d particles\n",p->n);
   for (int i=0; i<p->n; ++i) {
      const float npx = p->x[2*i+0] * (float)(pnx-1);
      const float npy = p->x[2*i+1] * (float)(pny-1) / yf;
      const float opx = p->oldx[2*i+0] * (float)(pnx-1);
      const float opy = p->oldx[2*i+1] * (float)(pny-1) / yf;

      // how many segments will we need?
      const int nseg = (int)(1.0 + 1.2 * sqrtf(powf(npx-opx,2)+powf(npy-opy,2)));
      //fprintf(stderr,"part %d needs %d segments\n", i, nseg); fflush(stderr);

      // drop the color on weighted by the ballistic coefficient (so heavy particles draw more heavily)
      //const float wgt = p->m[i];
      // weight by velocity also (04b)
      const float velmag = sqrt(pow(p->u[2*i+0],2)+pow(p->u[2*i+1],2)+1.e-6);
      //const float wgt = p->m[i] * velmag;
      // shimmer a little by reflecting off a normal direction, and account for vel (04c)
      //const float sheer = 10.0*(p->u[2*i+0]*0.8 + p->u[2*i+1]*0.2);
      //const float wgt = 0.5 * p->m[i] * (sheer+velmag);
      // shimmer more (04e)
      //const float sheer = (p->u[2*i+0]*0.95 - p->u[2*i+1]*0.1) / velmag;
      //const float wgt = (0.01*sheer*sheer*sheer*sheer + 0.5*velmag) / (float)nseg;
      const float wgt = draw_fac * powf(p->m[i], mass_pow) * powf(velmag, vel_pow) / (float)nseg;
      //if (i==0) fprintf(stdout,"  particle 1 is at %g %g with weight %g\n",npx,npy,wgt);

      for (int j=0; j<nseg; ++j) {
         const float tfac = (j+0.5) / (float)nseg;
         const float px = opx + tfac*(npx-opx);
         const float py = opy + tfac*(npy-opy);
         //fprintf(stderr,"  seg %d at %g %g\n", j, px, py); fflush(stderr);

         int ix = rintf(px-0.5);
         int iy = rintf(py-0.5);
         int ixp1 = ix+1;
         int iyp1 = iy+1;
         const float fx = px - (float)ix;
         const float fy = py - (float)iy;

         if (ix<0) ix += pnx-1;
         if (ixp1<0) ixp1 += pnx-1;
         if (iy<0) iy += pny-1;
         if (iyp1<0) iyp1 += pny-1;

         if (ix>=pnx) ix -= pnx-1;
         if (ixp1>=pnx) ixp1 -= pnx-1;
         if (iy>=pny) iy -= pny-1;
         if (iyp1>=pny) iyp1 -= pny-1;

         const float fac1 = wgt * (1.0-fx)*(1.0-fy);
         outred[ix][iy] = (outred[ix][iy] + fac1*p->c[3*i+0]) / (1.0+fac1);
         outgrn[ix][iy] = (outgrn[ix][iy] + fac1*p->c[3*i+1]) / (1.0+fac1);
         outblu[ix][iy] = (outblu[ix][iy] + fac1*p->c[3*i+2]) / (1.0+fac1);

         const float fac2 = wgt * fx*(1.0-fy);
         outred[ixp1][iy] = (outred[ixp1][iy] + fac2*p->c[3*i+0]) / (1.0+fac2);
         outgrn[ixp1][iy] = (outgrn[ixp1][iy] + fac2*p->c[3*i+1]) / (1.0+fac2);
         outblu[ixp1][iy] = (outblu[ixp1][iy] + fac2*p->c[3*i+2]) / (1.0+fac2);

         const float fac3 = wgt * (1.0-fx)*fy;
         outred[ix][iyp1] = (outred[ix][iyp1] + fac3*p->c[3*i+0]) / (1.0+fac3);
         outgrn[ix][iyp1] = (outgrn[ix][iyp1] + fac3*p->c[3*i+1]) / (1.0+fac3);
         outblu[ix][iyp1] = (outblu[ix][iyp1] + fac3*p->c[3*i+2]) / (1.0+fac3);

         const float fac4 = wgt * fx*fy;
         outred[ixp1][iyp1] = (outred[ixp1][iyp1] + fac4*p->c[3*i+0]) / (1.0+fac4);
         outgrn[ixp1][iyp1] = (outgrn[ixp1][iyp1] + fac4*p->c[3*i+1]) / (1.0+fac4);
         outblu[ixp1][iyp1] = (outblu[ixp1][iyp1] + fac4*p->c[3*i+2]) / (1.0+fac4);
      }
   }
}

//
// multiply the particle image by the mask image
//
// input mask image is nx by ny
// particle color image is pnx by pny
//
void mult_part_by_mask (float yf, int nx, int ny, float** mask,
                        int pnx, int pny, float** red, float** grn, float** blu) {

   float** tempin[1];
   tempin[0] = mask;
   float outvals[1];

   // copy the in colors to the out array, and scale both
   for (int ix=0; ix<pnx; ix++) {
      for (int iy=0; iy<pny; iy++) {
         const float px = ix / (float)(pnx-1);
         const float py = yf * iy / (float)(pny-1);
         // pick color from the input map at nx by ny (use the interp array call?)
         //const float maskval = interpolate_using_M4p_2d(nx,ny,PERIODIC,PERIODIC,mask,px,py);
         interpolate_array_using_TSC_2d(nx,ny,PERIODIC,PERIODIC,tempin,px,py,1,outvals);
         // apply it to the output map at pnx by pny
         red[ix][iy] *= outvals[0];
         grn[ix][iy] *= outvals[0];
         blu[ix][iy] *= outvals[0];
      }
   }
}

