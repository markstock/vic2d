/*
 * vic3d
 *
 * Copyright 2004,2005 Mark J. Stock mstock@umich.edu
 *
 * a 3D vortex method which uses the method of characteristics for the
 * convection step, and a single explicit step for diffusion and vorticity
 * creation from walls
 *
 * other files needed: gr3.c and f2c.h
 *
 * gcc -o vic3d -lf2c -lm vic3d.c libvicmoc.c gr23.c
 * gcc -o vic3d -O3 -march=pentium4 -mno-ieee-fp -ffast-math -lf2c -lm vic3d.c libvicmoc.c gr23.c
 * gcc -o vic3d -O3 -march=athlon-xp -mno-ieee-fp -ffast-math -lf2c -lm vic3d.c libvicmoc.c gr23.c
 * make
 */

#include "vicmoc.h"

int main(int argc,char **argv) {

   int isStam = FALSE;			// set to TRUE to run Stable Fluids calc
   int i,j,k,cnt;
   int nx,ny,nz;
   int xbdry,ybdry,zbdry;
   int ix,iy,iz;
   int step = 0;
   int outscale;
   int writeOutput = TRUE;
   int maxStep;

   int sc_cnt = 3;			// all 3 vorticities
   int sc_type[MAX_SCALARS];
   float sc_diffus[MAX_SCALARS];
   int use_SF = FALSE;
   int use_DF = FALSE;
   int sfi;

   // the field values
   float courant;
   float reynolds;
   float dt;
   float bn;
   float coeff,Re,yf,zf;
   float px,py,pz,temp,tottemp,avgtemp = 0.0;
   float ***u[3];			// all velocities, vorticity, etc
   float ***a[MAX_SCALARS];		// all velocities, vorticity, etc
   // float ***t[MAX_SCALARS];		// all velocities, vorticity, etc
   float ke = 0.;
   float ke_last = 0.;
   float enst = 0.;
   float nu_eff;

   // now for the particles
   int useParticles = FALSE;
   int pnum = 0;
   float **ploc;
   float **pvel;
   float *prad;

   // lastly, the bookkeeping
   unsigned long int tics,last_tics;	// timers
   char outfileroot[70];
   char progname[80];


   // set default simulation properties
   nx = 33;
   ny = 33;
   nz = 33;
   xbdry = WALL;
   ybdry = WALL;
   zbdry = WALL;
   courant = 10.0;
   reynolds = 1000.0;
   outscale = 1;
   bn = 1.;
   maxStep = 999;
   // only necessary scalar type is vorticity, always at position 0
   sc_type[0] = WX;
   sc_type[1] = WY;
   sc_type[2] = WZ;
   for (i=3; i<MAX_SCALARS; i++) sc_type[i] = -1;

   // read command-line
   (void) strcpy(progname,argv[0]);
   if (argc < 1) (void) Usage(progname,0);
   for (i=1; i<argc; i++) {
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
         use_SF = TRUE;
      } else if (strncmp(argv[i], "-re", 3) == 0) {
         reynolds = atof(argv[++i]);
      } else if (strncmp(argv[i], "-o", 2) == 0) {
         outscale = atoi(argv[++i]);
      } else if (strncmp(argv[i], "-cn", 3) == 0) {
         courant = atof(argv[++i]);
      } else if (strncmp(argv[i], "-b", 2) == 0) {
         bn = atof(argv[++i]);
      } else if (strncmp(argv[i], "-stam", 5) == 0) {
         isStam = TRUE;
      } else if (strncmp(argv[i], "-step", 5) == 0) {
         maxStep = atoi(argv[++i]);
      } else if (strncmp(argv[i], "-noout", 6) == 0) {
         writeOutput = FALSE;
      } else
         (void) Usage(progname,0);
   }

   // set y max bound (if square, yf=1.0, otherwise scale it as xmax=1.0
   yf = (float)(ny-1)/(float)(nx-1);
   zf = (float)(nz-1)/(float)(nx-1);

   // set time step (dt), base it on the maximum courant number (max speed of)
   dt = courant/(nx-1);

   // allocate arrays
   // three velocities
   u[XV] = allocate_3d_array_f(nx,ny,nz);
   u[YV] = allocate_3d_array_f(nx,ny,nz);
   u[ZV] = allocate_3d_array_f(nx,ny,nz);
   // and all scalars (3 vorticity plus others)
   a[WX] = allocate_3d_array_f(nx,ny,nz);
   a[WY] = allocate_3d_array_f(nx,ny,nz);
   a[WZ] = allocate_3d_array_f(nx,ny,nz);
   // t[WX] = allocate_3d_array_f(nx,ny,nz);
   // t[WY] = allocate_3d_array_f(nx,ny,nz);
   // t[WZ] = allocate_3d_array_f(nx,ny,nz);
   for (i=0; i<3; i++) sc_diffus[i] = 1.0/reynolds;
   if (use_SF) {
      sc_type[sc_cnt] = SF;
      sc_diffus[sc_cnt] = 1.0/reynolds;
      a[sc_cnt] = allocate_3d_array_f(nx,ny,nz);
      // t[sc_cnt] = allocate_3d_array_f(nx,ny,nz);
      sc_cnt++;
   }
   if (use_DF) {
      sc_type[sc_cnt] = DF;
      sc_diffus[sc_cnt] = 1.0/reynolds;
      a[sc_cnt] = allocate_3d_array_f(nx,ny,nz);
      // t[sc_cnt] = allocate_3d_array_f(nx,ny,nz);
      sc_cnt++;
   }

   // set initial vorticity (position 0)
   for (ix=0; ix<nx; ix++) {
      for (iy=0; iy<ny; iy++) {
         for (iz=0; iz<nz; iz++) {
            a[WX][ix][iy][iz] = 0.0;
            a[WY][ix][iy][iz] = 0.0;
            a[WZ][ix][iy][iz] = 0.0;
         }
      }
   }
   if (FALSE) {
      // create a wavy shear layer of positive vorticity (pseudo-2d)
      for (ix=0; ix<nx; ix++) {
         px = (float)ix/(float)(nx-1);
         for (iy=0; iy<ny; iy++) {
            py = (float)iy/(float)(nx-1);
            for (iz=0; iz<nz; iz++) {
               pz = (float)iz/(float)(nx-1);
               a[WY][ix][iy][iz] = 10.0 - 80.0*fabs((zf/2.0) - pz - 0.05*sin(px*2.0*M_PI));
               if (a[WY][ix][iy][iz] < 0.0) a[WY][ix][iy][iz] = 0.0;
            }
         }
      }
   } else if (FALSE) {
      // create a wavy shear layer of positive vorticity (fully 3d)
      for (ix=0; ix<nx; ix++) {
         px = (float)ix/(float)(nx-1);
         for (iy=0; iy<ny; iy++) {
            py = (float)iy/(float)(nx-1);
            for (iz=0; iz<nz; iz++) {
               pz = (float)iz/(float)(nx-1);
               a[WY][ix][iy][iz] = 10.0 - 80.0*fabs((zf/2.0) - pz - 0.05*sin(px*2.0*M_PI) - 0.025*sin(py*2.0*M_PI));
               if (a[WY][ix][iy][iz] < 0.0) a[WY][ix][iy][iz] = 0.0;
            }
         }
      }
   } else if (FALSE) {
      // create a random field of vorticity
      for (ix=0; ix<nx; ix++) {
         for (iy=0; iy<ny; iy++) {
            for (iz=0; iz<nz; iz++) {
               a[WX][ix][iy][iz] = 2.0*rand()/RAND_MAX - 1.0;
               a[WY][ix][iy][iz] = 2.0*rand()/RAND_MAX - 1.0;
               a[WZ][ix][iy][iz] = 2.0*rand()/RAND_MAX - 1.0;
            }
         }
      }
   } else if (FALSE) {
      // add a blob, similar to Hill's spherical vortex?
      // add_hills_spherical_vortex(nx,ny,nz,xbdry,ybdry,zbdry,a,0.6,0.2,0.3,0.1,10.0);
   }
   if (FALSE) {
      // create the atmospheric boundary layer - NEEDS WORK
      for (ix=0; ix<nx; ix++) {
         for (iy=0; iy<ny; iy++) {
            for (iz=0; iz<nz; iz++) {
               if ((float)(iz)/(0.2*(nz-1)) < 1.)
                  a[WY][ix][iy][iz] = 10.0*(0.5*cos((M_PI*iz)/(0.2*(nz-1)))+0.5);
            }
         }
      }
   }
   if (FALSE) {
      // create three blobs of vorticity
      //add_smooth_circular_blob_3d(nx,ny,nz,xbdry,ybdry,zbdry,0,a[WX],0.3,0.3,0.1,-10.0);
      //add_smooth_circular_blob_3d(nx,ny,nz,xbdry,ybdry,zbdry,1,a[WY],0.4,0.16,0.1,10.0);
      //add_smooth_circular_blob_3d(nx,ny,nz,xbdry,ybdry,zbdry,1,a[WY],0.42,0.5,0.08,20.0);
      //add_smooth_circular_blob_3d(nx,ny,nz,xbdry,ybdry,zbdry,1,a[WY],0.58,0.5,0.04,-20.0);
      add_smooth_circular_blob_3d(nx,ny,nz,xbdry,ybdry,zbdry,1,a[WY],0.44,0.5,0.08,20.0);
      add_smooth_circular_blob_3d(nx,ny,nz,xbdry,ybdry,zbdry,1,a[WY],0.56,0.5,0.04,-20.0);
   }
   if (TRUE) {
      // create a thick vortex ring
      // two oblique colliding rings:
      add_vortex_ring_3d(nx,ny,nz,xbdry,ybdry,zbdry,a[WX],a[WY],a[WZ],
         0.6,0.6,0.7,-1.,-1.,-1.,0.12,0.06,20.);
      add_vortex_ring_3d(nx,ny,nz,xbdry,ybdry,zbdry,a[WX],a[WY],a[WZ],
         0.4,0.4,0.7,1.,1.,-1.,0.12,0.06,20.);
      // two normal colliding rings:
      //add_vortex_ring_3d(nx,ny,nz,xbdry,ybdry,zbdry,a[WX],a[WY],a[WZ],
      //   0.2,0.5,0.5,1.,0.,0.,0.15,0.075,20.);
      //add_vortex_ring_3d(nx,ny,nz,xbdry,ybdry,zbdry,a[WX],a[WY],a[WZ],
      //   0.8,0.5,0.52,-1.,0.,0.,0.15,0.075,20.);
   }

   // set initial scalar
   for (i=1; i<MAX_SCALARS; i++) if (sc_type[i] == SF) sfi = i;
   if (use_SF) {
      for (ix=0; ix<nx; ix++) {
         for (iy=0; iy<ny; iy++) {
            for (iz=0; iz<nz; iz++) {
               a[sfi][ix][iy][iz] = 0.0;
            }
         }
      }
      if (TRUE) {
         // create a sharp blob of scalar
         add_smooth_spherical_blob(nx,ny,nz,xbdry,ybdry,zbdry,a[sfi],0.5,0.5,0.81,0.15,1.0);
         add_smooth_spherical_blob(nx,ny,nz,xbdry,ybdry,zbdry,a[sfi],0.51,0.5,0.2,0.15,-1.0);
      }
      if (FALSE) {
         // create a Gaussian blob of scalar
         add_singular_blob_3d(nx,ny,nz,xbdry,ybdry,zbdry,a[sfi],0.45,0.5,0.05,0.3,-0.00005);
      }
      if (FALSE) {
         // create a cube of scalar
         add_cube_3d(nx,ny,nz,xbdry,ybdry,zbdry,a[sfi],0.45,0.45,0.0,0.1,1.0);
      }
      if (FALSE) {
         // columns of scalar
         //for (ix=nx/2; ix<nx; ix++) {
         for (ix=0; ix<nx; ix++) {
            for (iy=0; iy<ny; iy++) {
               for (iz=0; iz<nz; iz++) {
                  //a[sfi][ix][iy][iz] = 1.0;
                  a[sfi][ix][iy][iz] = sin(ix*0.7)*sin(iy*0.5) +
                                       0.2*rand()/RAND_MAX - 0.1;
               }
            }
         }
      }
      if (FALSE) {
         // one column of scalar
         for (ix=0; ix<nx; ix++) {
            px = (float)ix/(float)(nx-1);
            for (iy=0; iy<ny; iy++) {
               py = (float)iy/(float)(ny-1);
               for (iz=0; iz<nz; iz++) {
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
      for (i=0;i<pnum;i++) ploc[i]=(float *)malloc(sizeof(float)*3);
      pvel = (float **)malloc(sizeof(float *)*pnum);
      for (i=0;i<pnum;i++) pvel[i]=(float *)malloc(sizeof(float)*3);
      prad = (float *)malloc(sizeof(float)*pnum);

      // place them
      cnt = 0;
      for (ix=0; ix<20; ix++) {
         for (iy=0; iy<20; iy++) {
            for (iz=0; iz<20; iz++) {
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

   // project the vorticity onto a divergence-free field
   make_solenoidal_3d(nx,ny,nz,xbdry,ybdry,zbdry,a[WX],a[WY],a[WZ]);

   // iterate through time
   step = 0;
   tottemp = 0.0;
   while (TRUE) {

      // fprintf(stdout,"\nBegin step %d\n",step);

      // ----------------------------
      // write output
      if (writeOutput) {
         sprintf(outfileroot,"wy_%04d",step);
         write_output_3d(outfileroot,nx,ny,nz,a[WY],-10.0,20.0,outscale,FALSE);
         if (use_SF) {
            sprintf(outfileroot,"outi_%04d",step);
            write_output_3d(outfileroot,nx,ny,nz,a[sfi],-0.1,0.2,outscale,TRUE);
            sprintf(outfileroot,"out_%04d",step);
            write_output_3d(outfileroot,nx,ny,nz,a[sfi],-0.5,1.0,outscale,FALSE);
         }
         if (step%10 == 0) {
         //sprintf(outfileroot,"vel_%04d",step);
         //write_3d_vtk(outfileroot,nx,ny,nz,u[XV],u[YV],u[ZV]);
         sprintf(outfileroot,"vort_%04d",step);
         write_3d_vtk(outfileroot,nx,ny,nz,a[WX],a[WY],a[WZ]);
         }
         // sprintf(outfileroot,"u_%04d",step);
         // write_output_3d(outfileroot,nx,ny,nz,u[XV],-0.5,1.0,outscale,FALSE);
         // sprintf(outfileroot,"v_%04d",step);
         // write_output_3d(outfileroot,nx,ny,nz,u[YV],-0.5,1.0,outscale,FALSE);
         // sprintf(outfileroot,"w_%04d",step);
         // write_output_3d(outfileroot,nx,ny,nz,u[ZV],-0.5,1.0,outscale,FALSE);
         if (useParticles) {
            sprintf(outfileroot,"part_%04d",step);
            write_output_particles_rad(outfileroot,pnum,ploc,pvel,prad);
         }
      }

      // accept input (only if interactive scheme)

      // set the timer
      last_tics = clock();

      // take one computational step forward in time
      Re = step_forward_3d(step,isStam,nx,ny,nz,xbdry,ybdry,zbdry,u,a,sc_cnt,sc_type,sc_diffus,dt,bn);

      // calculate energy and enstrophy
      ke_last = ke;
      ke = find_energy_3d(nx,ny,nz,xbdry,ybdry,zbdry,u[XV],u[YV],u[ZV]);
      enst = find_energy_3d(nx,ny,nz,xbdry,ybdry,zbdry,a[WX],a[WY],a[WZ]);
      printf("  ke %g  enst %g\n",ke,enst);

      // find effective diffusion
      nu_eff = (ke_last-ke)/(2.*enst*dt);
      printf("  effective nu %g  explicit nu %g\n",nu_eff,1./reynolds);

      // move the particles
      if (useParticles) {
        explicit_particle_move_3d(nx,ny,nz,xbdry,ybdry,zbdry,u,dt,Re,pnum,ploc,pvel);
      }

      // calculate the total time elapsed, and the time for this last step
      if (writeOutput) {
      tics = clock();
      temp = ((float)tics-(float)last_tics)/CLOCKS_PER_SEC;
      if (temp < 0.0) temp += 4.2949673e+09/(float)CLOCKS_PER_SEC;
      tottemp += temp;
      avgtemp = tottemp/(step+1);
      // fprintf(stdout,"  took %g sec\n",temp);
      // fprintf(stdout,"  total %g and mean %g sec\n",tottemp,avgtemp);
      fprintf(stdout,"mean step time: %.3g sec, effective Re %g\n",avgtemp,Re);
      }

      // ----------------------------
      // check end conditions
      step++;
      if (step > maxStep) break;
   }

   fprintf(stderr,"\nDone.\n");
   exit(0);
}


/*
 * This function writes basic usage information to stderr,
 * and then quits. Too bad.
 */
int Usage(char progname[80],int status) {

   static char **cpp, *help_message[] = {
   "where [-options] are one or more of the following:                         ",
   "                                                                           ",
   "   -x [int]    resolution in x-direction, default=32                       ",
   "                                                                           ",
   "   -y [int]    resolution in y-direction, default=32                       ",
   "                                                                           ",
   "   -z [int]    resolution in z-direction, default=32                       ",
   "                                                                           ",
   "   -px         make domain periodic in x-direction, default is not periodic",
   "                                                                           ",
   "   -py         make domain periodic in y-direction, default is not periodic",
   "                                                                           ",
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

