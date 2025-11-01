/*
 * vic3d.c - driver program for three-dimensional vortex-in-cell method of characteristics
 *
 * Copyright 2004,05,20,25 Mark J. Stock <markjstock@gmail.com>
 *
 * a 3D vortex method which uses the method of characteristics for the
 * convection step, and a single explicit step for diffusion and vorticity
 * creation from walls
 */

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <sys/time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "utility.h"
#include "inout.h"
#include "vicmoc.h"
#include "libvicmoc3d.h"

float compute_and_write_stats(const int, const int, const float, const float, const float, const double, const int, const int, const int, const int, const int, const int, float****, float****);
int Usage(char[MAXCHARS],int);

int main(int argc,char **argv) {

   //static bool first_time = true;
   bool isStam = false;
   int nx = 33;
   int ny = 33;
   int nz = 33;
   int xbdry,ybdry,zbdry;
   int step = 0;
   int ccnidx = 0;
   int outscale = 1;
   int writeOutput = true;
   int maxstep = 9999;
   int writeevery = 1;
   float vmax = 1.f;
   float simtime = 0.f;
   //bool use_16bpp = false;
   bool silent = true;

   int sc_cnt = 3;			// all 3 vorticities
   int sc_type[MAX_SCALARS];
   float sc_diffus[MAX_SCALARS];
   int use_SF = false;
   int use_DF = false;
   int sfi;

   // the field values
   float Re = 1000.f;
   float dt = -1.f;
   float bn = 1.f;
   float effRe = -1.f;
   //float yf,zf;
   float px,py,pz,temp;
   float courant = 10.f;
   float courantconst = -1.f;
   float *gravity;

   float ***u[3];			// all velocities, vorticity, etc
   float ***a[MAX_SCALARS];	// all velocities, vorticity, etc
   float ***t[MAX_SCALARS];	// all velocities, vorticity, etc
   float vortscale,velscale,scalescale;
   float ccnvmax[VMAXAVG];

   // now for the particles
   int useParticles = false;
   int pnum = 0;
   float **ploc;
   float **pvel;
   float *prad;

   // bookkeeping
   struct timeval t_curr, t_last;
   double walltime = 0.;
   float **tempArry;
   char outfileroot[MAXCHARS];
   char progname[MAXCHARS];


   // set default simulation properties
   xbdry = WALL;
   ybdry = WALL;
   zbdry = WALL;
   gravity = allocate_1d_array_f(3);
   gravity[0] = 0.0;
   gravity[1] = 0.0;
   gravity[2] = 1.0;
   vortscale = 100.0;
   velscale = 1.0;
   scalescale = 1.0;
   writeevery = 1;

   // only necessary scalar type is vorticity, always at position 0
   sc_type[0] = WX;
   sc_type[1] = WY;
   sc_type[2] = WZ;
   for (int i=3; i<MAX_SCALARS; i++) sc_type[i] = -1;

   // read command-line
   (void) strcpy(progname,argv[0]);
   if (argc < 1) (void) Usage(progname,0);
   for (int i=1; i<argc; i++) {

      if (strncmp(argv[i], "-x", 2) == 0) {
         nx = atoi(argv[++i]) + 1;
      } else if (strncmp(argv[i], "-y", 2) == 0) {
         ny = atoi(argv[++i]) + 1;
      } else if (strncmp(argv[i], "-z", 2) == 0) {
         nz = atoi(argv[++i]) + 1;
      } else if (strncmp(argv[i], "-px", 3) == 0) {
         xbdry = PERIODIC;
      } else if (strncmp(argv[i], "-py", 3) == 0) {
         ybdry = PERIODIC;
      } else if (strncmp(argv[i], "-pz", 3) == 0) {
         zbdry = PERIODIC;
      } else if (strncmp(argv[i], "-sc", 3) == 0) {
         use_SF = true;
      } else if (strncmp(argv[i], "-re", 3) == 0) {
         Re = atof(argv[++i]);
      } else if (strncmp(argv[i], "-o", 2) == 0) {
         outscale = atoi(argv[++i]);
      } else if (strncmp(argv[i], "-cn", 3) == 0) {
         courant = atof(argv[++i]);
      } else if (strncmp(argv[i], "-constcn", 8) == 0) {
         courantconst = atof(argv[++i]);
      } else if (strncmp(argv[i], "-b", 2) == 0) {
         bn = atof(argv[++i]);
      } else if (strncmp(argv[i], "-stam", 5) == 0) {
         isStam = true;
      } else if (strncmp(argv[i], "-dt", 3) == 0) {
         dt = atof(argv[++i]);
      } else if (strncmp(argv[i], "-every", 3) == 0) {
         writeevery = atoi(argv[++i]);
         if (writeevery < 1) writeevery = 1;
      } else if (strncmp(argv[i], "-step", 5) == 0) {
         maxstep = atoi(argv[++i]);
      } else if (strncmp(argv[i], "-noout", 6) == 0) {
         writeOutput = false;
      } else {
         fprintf(stderr,"Unknown option (%s)\n",argv[i]);
         (void) Usage(progname,0);
      }
   }

   // set y max bound (if square, yf=1.0, otherwise scale it as xmax=1.0
   //yf = (float)(ny-1)/(float)(nx-1);
   //zf = (float)(nz-1)/(float)(nx-1);

   // set time step (dt), base it on the maximum courant number (max speed of)
   //  but only if courant number is not set
   if (dt < 0.)
     dt = courant/(nx-1);
   else
     courant = dt*(nx-1);
   fprintf(stdout,"Using dt=%g, or courant number=%g\n",dt,courant);

   // -----------------------------------------------
   // Allocate and initialize

   // three velocities first
   u[XV] = allocate_3d_array_f(nx,ny,nz);
   u[YV] = allocate_3d_array_f(nx,ny,nz);
   u[ZV] = allocate_3d_array_f(nx,ny,nz);
   // and all scalars (3 vorticity plus others)
   a[WX] = allocate_3d_array_f(nx,ny,nz);
   a[WY] = allocate_3d_array_f(nx,ny,nz);
   a[WZ] = allocate_3d_array_f(nx,ny,nz);
   t[WX] = allocate_3d_array_f(nx,ny,nz);
   t[WY] = allocate_3d_array_f(nx,ny,nz);
   t[WZ] = allocate_3d_array_f(nx,ny,nz);
   // origin points from moc
   //xo[XV] = allocate_3d_array_f(nx,ny,nz);
   //xo[YV] = allocate_3d_array_f(nx,ny,nz);
   //xo[ZV] = allocate_3d_array_f(nx,ny,nz);
   for (int i=0; i<3; i++) sc_diffus[i] = 1.0/Re;
   if (use_SF) {
      sc_type[sc_cnt] = SF;
      sc_diffus[sc_cnt] = 1.0/Re;
      a[sc_cnt] = allocate_3d_array_f(nx,ny,nz);
      t[sc_cnt] = allocate_3d_array_f(nx,ny,nz);
      sc_cnt++;
   }
   if (use_DF) {
      sc_type[sc_cnt] = DF;
      sc_diffus[sc_cnt] = 1.0/Re;
      a[sc_cnt] = allocate_3d_array_f(nx,ny,nz);
      t[sc_cnt] = allocate_3d_array_f(nx,ny,nz);
      sc_cnt++;
   }

   // allocate space for temporary arrays
   //for (int i=0; i<MAX_SCALARS; i++) {
   //   if (a[i] != NULL) {
   //      t[i] = allocate_2d_array_f(nx,ny);
   //      if (!silent) fprintf(stdout,"Array %d has diffusivity %g\n",i,sc_diffus[i]);
   //   }
   //}

   // -----------------------------------------------
   // Set initial conditions

   // Set vorticity ----------------------------------

   // zero initial vorticity
   for (int ix=0; ix<nx; ix++) {
      for (int iy=0; iy<ny; iy++) {
         for (int iz=0; iz<nz; iz++) {
            a[WX][ix][iy][iz] = 0.f;
            a[WY][ix][iy][iz] = 0.f;
            a[WZ][ix][iy][iz] = 0.f;
         }
      }
   }
   if (false) {
      // create a wavy shear layer of positive vorticity (pseudo-2d)
      const float zf = (float)(nz-1)/(float)(nx-1);
      for (int ix=0; ix<nx; ix++) {
         px = (float)ix/(float)(nx-1);
         for (int iy=0; iy<ny; iy++) {
            py = (float)iy/(float)(nx-1);
            for (int iz=0; iz<nz; iz++) {
               pz = (float)iz/(float)(nx-1);
               a[WY][ix][iy][iz] = 10.0 - 80.0*fabs((zf/2.0) - pz - 0.05*sin(px*2.0*M_PI));
               if (a[WY][ix][iy][iz] < 0.0) a[WY][ix][iy][iz] = 0.0;
            }
         }
      }
   } else if (false) {
      // create a wavy shear layer of positive vorticity (fully 3d)
      const float zf = (float)(nz-1)/(float)(nx-1);
      for (int ix=0; ix<nx; ix++) {
         px = (float)ix/(float)(nx-1);
         for (int iy=0; iy<ny; iy++) {
            py = (float)iy/(float)(nx-1);
            for (int iz=0; iz<nz; iz++) {
               pz = (float)iz/(float)(nx-1);
               a[WY][ix][iy][iz] = 10.0 - 80.0*fabs((zf/2.0) - pz - 0.05*sin(px*2.0*M_PI) - 0.025*sin(py*2.0*M_PI));
               if (a[WY][ix][iy][iz] < 0.0) a[WY][ix][iy][iz] = 0.0;
            }
         }
      }
   } else if (false) {
      // create a random field of vorticity
      for (int ix=0; ix<nx; ix++) {
         for (int iy=0; iy<ny; iy++) {
            for (int iz=0; iz<nz; iz++) {
               a[WX][ix][iy][iz] = 2.0*rand()/RAND_MAX - 1.0;
               a[WY][ix][iy][iz] = 2.0*rand()/RAND_MAX - 1.0;
               a[WZ][ix][iy][iz] = 2.0*rand()/RAND_MAX - 1.0;
            }
         }
      }
   } else if (false) {
      // add a blob, similar to Hill's spherical vortex?
      // add_hills_spherical_vortex(nx,ny,nz,xbdry,ybdry,zbdry,a,0.6,0.2,0.3,0.1,10.0);
   }

   if (false) {
      // create a tube of vorticity along x
      for (int ix=0; ix<nx; ix++) {
         for (int iy=0; iy<ny; iy++) {
            py = (float)iy/(float)(nx-1);
            for (int iz=0; iz<nz; iz++) {
               pz = (float)iz/(float)(nx-1);
               px = sqrt(pow(py-0.5,2) + pow(pz-0.3,2));
               if (px < 0.1) {
                  a[WX][ix][iy][iz] = 1.0;
               } else if (px < 0.2) {
                  a[WX][ix][iy][iz] = 0.5*(cos(M_PI*(px-0.1)/0.1)+1.);
               }
            }
         }
      }
   }
   if (false) {
      // create a tube of vorticity along x
      for (int ix=0; ix<nx; ix++) {
         for (int iy=0; iy<ny; iy++) {
            py = (float)iy/(float)(nx-1);
            for (int iz=0; iz<nz; iz++) {
               pz = (float)iz/(float)(nx-1);
               px = sqrt(pow(py-0.4,2) + pow(pz-0.7,2));
               if (px < 0.1) {
                  a[WX][ix][iy][iz] = -1.0;
               } else if (px < 0.2) {
                  a[WX][ix][iy][iz] = -0.5*(cos(M_PI*(px-0.1)/0.1)+1.);
               }
            }
         }
      }
   }
   if (false) {
      // create the atmospheric boundary layer - NEEDS WORK
      for (int ix=0; ix<nx; ix++) {
         for (int iy=0; iy<ny; iy++) {
            for (int iz=0; iz<nz; iz++) {
               if ((float)(iz)/(0.2*(nz-1)) < 1.)
                  a[WY][ix][iy][iz] = 10.0*(0.5*cos((M_PI*iz)/(0.2*(nz-1)))+0.5);
            }
         }
      }
   }
   if (false) {
      // create three blobs of vorticity
      //add_smooth_circular_blob_3d(nx,ny,nz,xbdry,ybdry,zbdry,0,a[WX],0.3,0.3,0.1,-10.0);
      //add_smooth_circular_blob_3d(nx,ny,nz,xbdry,ybdry,zbdry,1,a[WY],0.4,0.16,0.1,10.0);
      //add_smooth_circular_blob_3d(nx,ny,nz,xbdry,ybdry,zbdry,1,a[WY],0.42,0.5,0.08,20.0);
      //add_smooth_circular_blob_3d(nx,ny,nz,xbdry,ybdry,zbdry,1,a[WY],0.58,0.5,0.04,-20.0);
      add_smooth_circular_blob_3d(nx,ny,nz,xbdry,ybdry,zbdry,1,a[WY],0.44,0.5,0.08,20.0);
      add_smooth_circular_blob_3d(nx,ny,nz,xbdry,ybdry,zbdry,1,a[WY],0.56,0.5,0.04,-20.0);
   }
   if (false) {
      // two oblique colliding rings:
      add_vortex_ring_3d(nx,ny,nz,xbdry,ybdry,zbdry,a[WX],a[WY],a[WZ],0.6,0.6,0.7,-1.,-1.,-1.,0.12,0.06,20.);
      add_vortex_ring_3d(nx,ny,nz,xbdry,ybdry,zbdry,a[WX],a[WY],a[WZ],0.4,0.4,0.7,1.,1.,-1.,0.12,0.06,20.);
   }
   if (true) {
      // two normal colliding rings:
      //add_vortex_ring_3d(nx,ny,nz,xbdry,ybdry,zbdry,a[WX],a[WY],a[WZ],0.4,0.5,0.5,1.,0.,0.,0.15,0.075,20.);
      //add_vortex_ring_3d(nx,ny,nz,xbdry,ybdry,zbdry,a[WX],a[WY],a[WZ],0.6,0.5,0.52,-1.,0.,0.,0.15,0.075,20.);
      //add_vortex_ring_3d(nx,ny,nz,xbdry,ybdry,zbdry,a[WX],a[WY],a[WZ],0.5,0.5,0.4,0.f,0.f,1.f,0.2,0.075,177.);
      //add_vortex_ring_3d(nx,ny,nz,xbdry,ybdry,zbdry,a[WX],a[WY],a[WZ],0.5,0.5,0.6,0.f,0.f,-1.f,0.2,0.075,177.);
      add_vortex_ring_3d(nx,ny,nz,xbdry,ybdry,zbdry,a[WX],a[WY],a[WZ],0.5,0.62,0.4,0.f,0.f,1.f,0.2,0.075,177.);
      add_vortex_ring_3d(nx,ny,nz,xbdry,ybdry,zbdry,a[WX],a[WY],a[WZ],0.5,0.58,0.6,0.f,0.f,-1.f,0.2,0.075,177.);
   }

   // set initial scalar
   for (int i=1; i<MAX_SCALARS; i++) if (sc_type[i] == SF) sfi = i;
   if (use_SF) {
      for (int ix=0; ix<nx; ix++) {
         for (int iy=0; iy<ny; iy++) {
            for (int iz=0; iz<nz; iz++) {
               a[sfi][ix][iy][iz] = 0.0;
            }
         }
      }
      if (false) {
         // create a sharp blob of scalar
         add_smooth_spherical_blob(nx,ny,nz,xbdry,ybdry,zbdry,a[sfi],0.5,0.5,0.81,0.15,1.0);
         add_smooth_spherical_blob(nx,ny,nz,xbdry,ybdry,zbdry,a[sfi],0.51,0.5,0.2,0.15,-1.0);
      }
      if (false) {
         // create a Gaussian blob of scalar
         add_singular_blob_3d(nx,ny,nz,xbdry,ybdry,zbdry,a[sfi],0.45,0.5,0.05,0.3,-0.00005);
      }
      if (false) {
         // create a cube of scalar
         add_cube_3d(nx,ny,nz,xbdry,ybdry,zbdry,a[sfi],0.45,0.45,0.0,0.1,1.0);
      }
      if (true) {
         // columns of scalar
         for (int ix=0; ix<nx; ix++) {
            for (int iy=0; iy<ny; iy++) {
               for (int iz=0; iz<nz; iz++) {
                  pz = (float)iz/(float)(nz-1);
                  a[sfi][ix][iy][iz] = 4.8*MAX(0.,-cos(2.*M_PI*pz) + 0.01*rand()/RAND_MAX - 0.8);
               }
            }
         }
      }
      if (false) {
         // columns of scalar
         //for (int ix=nx/2; ix<nx; ix++) {
         for (int ix=0; ix<nx; ix++) {
            for (int iy=0; iy<ny; iy++) {
               for (int iz=0; iz<nz; iz++) {
                  //a[sfi][ix][iy][iz] = 1.0;
                  a[sfi][ix][iy][iz] = sin(ix*0.7)*sin(iy*0.5) +
                                       0.2*rand()/RAND_MAX - 0.1;
               }
            }
         }
      }
      if (false) {
         // one column of scalar
         for (int ix=0; ix<nx; ix++) {
            px = (float)ix/(float)(nx-1);
            for (int iy=0; iy<ny; iy++) {
               py = (float)iy/(float)(ny-1);
               for (int iz=0; iz<nz; iz++) {
                  temp = sqrt(pow(0.5-px,2)+pow(0.5-py,2));
                  temp += 0.02*rand()/RAND_MAX - 0.01;
                  if (temp < 0.11) {
                     a[sfi][ix][iy][iz] = 0.5*(1.+cos((0.1-temp)*100.*M_PI));
                  }
                  if (temp < 0.09) {
                     a[sfi][ix][iy][iz] = 1.;
                  }
               }
            }
         }
      }
   }

   // set particles
   if (useParticles) {
      pnum = 20*20*20;
      // allocate memory for particles
      ploc = (float **)malloc(sizeof(float *)*pnum);
      for (int i=0;i<pnum;i++) ploc[i]=(float *)malloc(sizeof(float)*3);
      pvel = (float **)malloc(sizeof(float *)*pnum);
      for (int i=0;i<pnum;i++) pvel[i]=(float *)malloc(sizeof(float)*3);
      prad = (float *)malloc(sizeof(float)*pnum);

      // place them
      int cnt = 0;
      for (int ix=0; ix<20; ix++) {
         for (int iy=0; iy<20; iy++) {
            for (int iz=0; iz<20; iz++) {
               ploc[cnt][0] = 0.45+0.1*((float)(ix)+0.5)/20.;
               ploc[cnt][1] = 0.45+0.1*((float)(iy)+0.5)/20.;
               ploc[cnt][2] = 0.0+0.1*((float)(iz)+0.5)/20.;
               pvel[cnt][0] = 0.;
               pvel[cnt][1] = 0.;
               pvel[cnt][2] = 0.;
               prad[cnt] = 0.1/20.;
               cnt++;
            }
         }
      }
   }


   // Initial velocity solve -------------------------

   fprintf(stdout,"Initial velocity solution\n");

   // project the vorticity onto a divergence-free field
   //make_solenoidal_3d(nx,ny,nz,xbdry,ybdry,zbdry,a[WX],a[WY],a[WZ]);

   // if you want to show velocity on the first step, solve for it here
   (void) find_vels_3d(step,nx,ny,nz,xbdry,ybdry,zbdry,u[XV],u[YV],u[ZV],a[WX],a[WY],a[WZ]);

   // find vmax (might need this)
   vmax = find_vmax(u[XV],u[YV],u[ZV],nx,ny,nz);
   // initialize the vmax averaging array
   vmax = courantconst / (dt*(nx+1));
   for (int i=0; i<VMAXAVG; i++) ccnvmax[i] = vmax;


   // ----------------------------------

   // iterate through time
   step = 0;
   while (true) {

      fprintf(stdout,"\nBegin step %d\n",step);

      // ----------------------------
      // write output
      if (writeOutput && step%writeevery == 0) {
         if (true) {
            const float magscale = 0.2f * vortscale;
            sprintf(outfileroot,"wm_%04d",step);
            write_vortmag_3d(outfileroot,nx,ny,nz,a[WX],a[WY],a[WZ],0.f,magscale,outscale);
         }
         if (true) {
            const int midplane = nx/2;
            sprintf(outfileroot,"wx_%04d",step);
            //write_output_3d(outfileroot,nx,ny,nz,a[WY],-10.0,20.0,outscale,false);
            write_png (outfileroot,ny,nz,false,false,
                       a[WX][midplane],-vortscale,2.*vortscale,
                       NULL,0.0,1.0,
                       NULL,0.0,1.0);
         }
         if (false) {
            const int midplane = nx/2;
            sprintf(outfileroot,"wy_%04d",step);
            write_png (outfileroot,ny,nz,false,false,
                       a[WY][midplane],-vortscale,2.*vortscale,
                       NULL,0.0,1.0,
                       NULL,0.0,1.0);
            sprintf(outfileroot,"wz_%04d",step);
            write_png (outfileroot,ny,nz,false,false,
                       a[WZ][midplane],-vortscale,2.*vortscale,
                       NULL,0.0,1.0,
                       NULL,0.0,1.0);
         }

         if (false) {
            const int midplane = nx/2;
            sprintf(outfileroot,"ux_%04d",step);
            write_png (outfileroot,ny,nz,false,false,
                       u[XV][midplane],-velscale,2.*velscale,
                       NULL,0.0,1.0,
                       NULL,0.0,1.0);
         }
         if (true) {
            const int midplane = nx/2;
            sprintf(outfileroot,"uy_%04d",step);
            write_png (outfileroot,ny,nz,false,false,
                       u[YV][midplane],-velscale,2.*velscale,
                       NULL,0.0,1.0,
                       NULL,0.0,1.0);
            sprintf(outfileroot,"uz_%04d",step);
            write_png (outfileroot,ny,nz,false,false,
                       u[ZV][midplane],-velscale,2.*velscale,
                       NULL,0.0,1.0,
                       NULL,0.0,1.0);
         }
         if (use_SF) {
            //sprintf(outfileroot,"outi_%04d",step);
            //write_output_3d(outfileroot,nx,ny,nz,a[sfi],-0.1,0.2,outscale,true);
            //write_png (outfileroot,ny,nz,false,false,
            //           a[sfi],-scalescale,2.*scalescale,
            //           NULL,0.0,1.0,
            //           NULL,0.0,1.0);
            sprintf(outfileroot,"out_%04d",step);
            //write_output_3d(outfileroot,nx,ny,nz,a[sfi],-0.5,1.0,outscale,false);
            write_png (outfileroot,ny,nz,false,false,
                       a[sfi][3],-scalescale,2.*scalescale,
                       NULL,0.0,1.0,
                       NULL,0.0,1.0);
            // integrated image
            sprintf(outfileroot,"outi_%04d",step);
            tempArry = flatten_to_2d(a[sfi],0,nx,ny,nz);
            write_png (outfileroot,ny,nz,false,false,
                       tempArry,-scalescale,2.*scalescale,
                       NULL,0.0,1.0,
                       NULL,0.0,1.0);
            free_2d_array_f(tempArry);
         }
         //if (step%10 == 0) {
         if (false) {
            //sprintf(outfileroot,"vel_%04d",step);
            //write_3d_vtk(outfileroot,nx,ny,nz,u[XV],u[YV],u[ZV]);
            sprintf(outfileroot,"vort_%04d",step);
            write_3d_vtk(outfileroot,nx,ny,nz,a[WX],a[WY],a[WZ]);
         }
         // sprintf(outfileroot,"u_%04d",step);
         // write_output_3d(outfileroot,nx,ny,nz,u[XV],-0.5,1.0,outscale,false);
         // sprintf(outfileroot,"v_%04d",step);
         // write_output_3d(outfileroot,nx,ny,nz,u[YV],-0.5,1.0,outscale,false);
         // sprintf(outfileroot,"w_%04d",step);
         // write_output_3d(outfileroot,nx,ny,nz,u[ZV],-0.5,1.0,outscale,false);
         if (useParticles) {
            sprintf(outfileroot,"part_%04d",step);
            write_output_particles_rad(outfileroot,pnum,ploc,pvel,prad);
         }
      }

      // check end conditions
      if (step >= maxstep) break;

      // accept input (only if interactive scheme)

      // if we want a constant cn time step, compute it here
      if (courantconst > 0.0) {
         // average vmax over a number of steps
         ccnidx = (ccnidx+1)%VMAXAVG;
         ccnvmax[ccnidx] = vmax;
         // if any were zero, reset them
         for (int i=0; i<VMAXAVG; i++) if (ccnvmax[i] < 1.e-6) ccnvmax[i] = vmax;
         // find new thing
         vmax = 0.0;
         for (int i=0; i<VMAXAVG; i++) vmax += ccnvmax[i];
         vmax /= (float)VMAXAVG;
         // recalculate dt
         if ((courantconst/vmax)/(nx+1) < 1.5*dt) {
            dt = (courantconst/vmax)/(nx+1);
            fprintf(stderr,"  changing dt to %g\n",dt);
         } else {
            dt *= 1.5;
         }
      }

      // set the timer
      gettimeofday(&t_last, 0);

      // take one computational step forward in time
      effRe = step_forward_3d(step,isStam,
              nx,ny,nz,xbdry,ybdry,zbdry,
              u,a,t,
              sc_cnt,sc_type,sc_diffus,
              dt,
              bn);

      // move the particles
      if (useParticles) {
        explicit_particle_move_3d(nx,ny,nz,xbdry,ybdry,zbdry,u,dt,effRe,pnum,ploc,pvel);
      }

      // calculate the total time elapsed, and the time for this last step
      gettimeofday(&t_curr, 0);
      walltime = (t_curr.tv_sec - t_last.tv_sec) + (t_curr.tv_usec - t_last.tv_usec)*1e-6;

      //if (writeOutput) {
         //tottemp += walltime;
         //float avgtemp = tottemp/(step+1);
         // fprintf(stdout,"  took %g sec\n",temp);
         // fprintf(stdout,"  total %g and mean %g sec\n",tottemp,avgtemp);
         //fprintf(stdout,"  mean step time: %.3g sec, effective Re %g\n",avgtemp,effRe);
         //fprintf(stdout,"  step time: %.4g sec, effective Re %g\n", walltime, effRe);
      //}

      // write statistics
      if (writeOutput) {
         vmax = compute_and_write_stats(silent,step,dt,Re,simtime,walltime,nx,ny,nz,xbdry,ybdry,zbdry,u,a);
      }

      // ----------------------------
      step++;
      simtime += dt;
   }

   fprintf(stderr,"\nDone.\n");
   exit(0);
}


/*
 * Do just what it says
 */
float compute_and_write_stats(const int silent, const int step,
      const float dt, const float Re, const float simtime, const double walltime,
      const int nx, const int ny, const int nz, const int xbdry, const int ybdry, const int zbdry,
      float ****u, float ****a) {

   static float ke_last = -1.f;
   static double total_time = 0.0;
   static char outfile[MAXCHARS];
   static bool initialized = false;
   static FILE *outp;

   // if last time through, close file pointer and RETURN
   if (step == -1) {
      fclose(outp);
      return(0.f);
   }

   // if first time through, prep stuff
   if (!initialized) {
      sprintf(outfile,"stat.dat");
      outp = fopen(outfile,"w");
      if (outp == NULL) {
         fprintf(stderr,"Could not open input file %s\n",outfile);
         fflush(stderr);
         exit(0);
      }
      initialized = true;
      fprintf(outp,"# step, sim time, vel max, vort max, KE, enst, CNw, cpu time step, cpu time total\n");
      fflush(outp);
   }

   // find vmax
   const float vmax = find_vmax(u[XV],u[YV],u[ZV],nx,ny,nz);

   // find vort max
   const float wmax = find_vmax(a[WX],a[WY],a[WZ],nx,ny,nz);

   // compute CN
   const float cn = vmax*dt*(float)nx;
   const float cnw = wmax*dt*(float)nx;
   printf("  CNv %g  CNw %g\n", cn, cnw);

   // calculate energy and enstrophy
   const float ke = find_energy_3d(nx,ny,nz,xbdry,ybdry,zbdry,u[XV],u[YV],u[ZV]);
   const float enst = find_energy_3d(nx,ny,nz,xbdry,ybdry,zbdry,a[WX],a[WY],a[WZ]);
   //printf("  ke %g  enst %g\n",ke,enst);

   // find effective diffusion
   const float nu_eff = (ke_last < 0.f) ? 0.f : (ke_last-ke)/(2.*enst*dt);
   ke_last = ke;
   printf("  effective nu %g  explicit nu %g\n",nu_eff,1./Re);

   // and time
   total_time += walltime;

   // write the stats here
   if (!silent) {
      if (step%10 == 0) fprintf(stdout,"# step, sim time, vel max, vort max, KE, enst, CNw, cpu step time, cpu total time\n");
      fprintf(stdout,"%d %g %g %g %g %g %g %g %g\n",step,simtime,vmax,wmax,ke,enst,cnw,walltime,total_time);
      fflush(stdout);
   }
   fprintf(outp,"%d %g %g %g %g %g %g %g %g\n",step,simtime,vmax,wmax,ke,enst,cnw,walltime,total_time);
   fflush(outp);

   return(vmax);
}


/*
 * This function writes basic usage information to stderr,
 * and then quits. Too bad.
 */
int Usage(char progname[MAXCHARS],int status) {

   static char **cpp, *help_message[] = {
   "where [-options] are one or more of the following:                         ",
   "                                                                           ",
   "   -x [int]    resolution in x-direction, default=32                       ",
   "   -y [int]    resolution in y-direction, default=32                       ",
   "   -z [int]    resolution in z-direction, default=32                       ",
   "                                                                           ",
   "   -px         make domain periodic in x-direction, default is not periodic",
   "   -py         make domain periodic in y-direction, default is not periodic",
   "   -pz         make domain periodic in z-direction, default is not periodic",
   "                                                                           ",
   "   -sc         track scalar fraction field                                 ",
   "                                                                           ",
   "   -o num      set output resolution to num times simulation resolution,   ",
   "               default=1                                                   ",
   "                                                                           ",
   "   -re [float] Reynolds number for flow, default=1000                      ",
   "                                                                           ",
   "   -cn [float] set the Courant number, based on a speed of 1.0; try        ",
   "               10 to 50 for a small simulation (32 to 128) or 50 to 500    ",
   "               for a larger one (256 to 2048); default=10                  ",
   "                                                                           ",
   "   -dt [float] time step size (if not set, use cn instead)                 ",
   "                                                                           ",
   "   -b [float]  set the Boussinesq number; default=1                        ",
   "                                                                           ",
   "   -stam       use Stam's velocity method-of-characteristics               ",
   "                                                                           ",
   "   -help       (in place of infile) returns this help information          ",
   " ",
   "Options may be abbreviated to an unambiguous length.",
   "Output is to a series of PNG files.",
   NULL
   };

   fprintf(stderr, "usage:\n  %s [options]\n\n", progname);
   for (cpp = help_message; *cpp; cpp++) fprintf(stderr, "%s\n", *cpp);
   fflush(stderr);
   exit(status);
   return(0);
}

