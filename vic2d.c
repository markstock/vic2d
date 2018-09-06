/*
 * vic2d
 *
 * Copyright 2004-18 Mark J. Stock mstock@umich.edu
 *
 * a 2D vortex method which uses the method of characteristics for the
 * convection step, and a single explicit step for diffusion and vorticity
 * creation from walls
 *
 */

#include "vicmoc.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include <ctype.h>

float compute_and_write_stats(int, int, float, float, float, int, int, int, int, float***);
void paint_splat (float,float,float,float,float,float,float,float,float,int,int,int,float**,float**,float**);
int get_random_color (float***,int,int,float*);
int Usage(char[80],int);

int main(int argc,char **argv) {

   static int first_time = TRUE;
   int isStam;
   int i,j;
   int nx,ny;
   int cnx,cny;
   int xbdry,ybdry;
   int ix,iy;
   int step = -1;
   int maxstep = 0;
   int writeevery;
   int outscale;
   int use_16bpp = FALSE;
   int silent = FALSE;

   float sc_diffus[MAX_SCALARS];
   int use_TEMP = FALSE;
   int use_DF = FALSE;
   int use_COLOR = FALSE;
   int reread_color = FALSE;
   int reread_vort = FALSE;
   int darkenonly = TRUE;
   int freeze_flow = FALSE;
   int freeze_at = -1;
   int stop_flow = FALSE;
   int stop_at = -1;
   int use_MASK = FALSE;
   int use_image = FALSE;
   int use_strong_strat = FALSE;
   int overlay_color_lr = FALSE;
   int overlay_color_tb = FALSE;

   int writeOutput = TRUE;
   int print_vort = FALSE;
   int print_vel = FALSE;
   int print_temp = FALSE;
   int print_mu = FALSE;
   int print_mask = FALSE;

   float md,mdlow,mdhigh;		// dimensional momentum diffusivity
   float vd;				// dimensional viscosity diffusivity
   float td,tdlow,tdhigh;		// dimensional density diffusivity
   float cd,cdlow,cdhigh;		// dimensional color diffusivity
   float bn;				// Boussinesq coefficient
   float dens_ratio;			// ratio of max to min density
   float heat_coeff;			// coefficient on the heat term
   float dt;
   float courant,courantconst;		// non-dim time step size
   float *freestream;
   float *wallvel;
   int recalc_vel;
   int move_colors;
   int gravtype;
   int ccnidx = 0;
   float *gravity;
   float coeff,yf,thickness;
   float effective_re;
   float px,py,cx,cy,rad,temp,temp2,maxval,minval,scale;
   float vortscale,velscale,vmax,velsq;
   float ccnvmax[VMAXAVG];
   float randvortscale = -1;
   float overlay_fraction = 0.01;
   float cputime = 0.;
   float totcputime = 0.;
   float avgcputime = 0.;
   float simtime = 0.;
   float thisc[3];
   float maskerr;
   float **color_left, **color_right, **color_top, **color_bottom;

   float **u[2];			// velocities
   float **a[MAX_SCALARS];		// other scalars
   float **t[MAX_SCALARS];		// temporary velocities, vorticity, etc
   float **mask;			// flow mask
   float **heat;			// constant heat source map
   float **c[3];			// color image storage
   float **acc[2];			// Lagrangian acceleration
   float **shear;			// shear magnitude

   // bookkeeping
   unsigned long int tics,last_tics;
   char progname[160];
   char outfileroot[150];
   int use_vort_img = FALSE;
   char vortfilename[160];
   int use_temp_img = FALSE;
   char tempfilename[160];
   int use_heat_img = FALSE;
   char heatfilename[160];
   int use_color_img = FALSE;
   char colorfilename[160];
   int use_div_img = FALSE;
   char divfilename[160];
   int use_mask_img = FALSE;
   char maskfilename[160];
   char colorsrcfilename[160];
   char mdfilename[160];
   char tdfilename[160];
   char cdfilename[160];


   // set default simulation properties
   nx = 513;
   ny = 513;
   xbdry = WALL;
   ybdry = WALL;
   isStam = FALSE;
   courant = 10.0;
   courantconst = -1.0;
   dt = -1.0;
   md = 1.e-3;
   mdlow = -1.;
   mdhigh = -1.;
   vd = 1.e-3;
   cd = 1.e-3;
   cdlow = -1.;
   cdhigh = -1.;
   heat_coeff = 1.;
   outscale = 1;
   writeevery = 1;
   maxstep = 99999;
   bn = 1.;
   dens_ratio = 1.;
   for (i=0; i<MAX_SCALARS; i++) sc_diffus[i] = 1.e-3;
   for (i=0; i<MAX_SCALARS; i++) a[i] = NULL;
   for (i=0; i<MAX_SCALARS; i++) t[i] = NULL;
   writeOutput = TRUE;
   freestream = allocate_1d_array_f(2);
   freestream[0] = 0.0;
   freestream[1] = 0.0;
   wallvel = allocate_1d_array_f(4);
   wallvel[0] = 0.0;	// left
   wallvel[1] = 0.0;	// right
   wallvel[2] = 0.0;	// bottom
   wallvel[3] = 0.0;	// top
   recalc_vel = TRUE;
   move_colors = TRUE;
   gravtype = 0;
   gravity = allocate_1d_array_f(2);
   gravity[0] = 0.0;
   gravity[1] = 1.0;
   vortscale = 10.;
   velscale = 1.;
   maskerr = 1.e-3;

   // read command-line
   (void) strcpy(progname,argv[0]);
   if (argc < 1) (void) Usage(progname,0);
   for (i=1; i<argc; i++) {

      if (strncmp(argv[i], "-x", 2) == 0) {
         nx = atoi(argv[++i]) + 1;
      } else if (strncmp(argv[i], "-y", 2) == 0) {
         ny = atoi(argv[++i]) + 1;
      } else if (strncmp(argv[i], "-px", 3) == 0) {
         xbdry = PERIODIC;
      } else if (strncmp(argv[i], "-py", 3) == 0) {
         ybdry = PERIODIC;
      } else if (strncmp(argv[i], "-open", 3) == 0) {
         xbdry = OPEN;
         ybdry = OPEN;
      } else if (strncmp(argv[i], "-bcl", 4) == 0) {
         wallvel[0] = atof(argv[++i]);
      } else if (strncmp(argv[i], "-bcr", 4) == 0) {
         wallvel[1] = atof(argv[++i]);
      } else if (strncmp(argv[i], "-bcb", 4) == 0) {
         wallvel[2] = atof(argv[++i]);
      } else if (strncmp(argv[i], "-bct", 4) == 0) {
         wallvel[3] = atof(argv[++i]);
      } else if (strncmp(argv[i], "-cn", 3) == 0) {
         courant = atof(argv[++i]);

      } else if (strncmp(argv[i], "-vf", 3) == 0) {
         strcpy (vortfilename,argv[++i]);
         use_vort_img = TRUE;
      } else if (strncmp(argv[i], "-vr", 3) == 0) {
         reread_vort = TRUE;
      } else if (strncmp(argv[i], "-vscale", 3) == 0) {
         vortscale = atof(argv[++i]);
      } else if (strncmp(argv[i], "-vdf", 4) == 0) {
         strcpy (mdfilename,argv[++i]);
         mdlow = atof(argv[++i]);
         mdhigh = atof(argv[++i]);
         vd = atof(argv[++i]);
      } else if (strncmp(argv[i], "-muprint", 3) == 0) {
         print_mu = TRUE;
      } else if (strncmp(argv[i], "-vd", 3) == 0) {
         md = atof(argv[++i]);
      } else if (strncmp(argv[i], "-vprint", 3) == 0) {
         print_vort = TRUE;
      } else if (strncmp(argv[i], "-uscale", 3) == 0) {
         velscale = atof(argv[++i]);
      } else if (strncmp(argv[i], "-uprint", 3) == 0) {
         print_vel = TRUE;
      } else if (strncmp(argv[i], "-freeze", 3) == 0) {
         freeze_flow = TRUE;
         if (argc > i+1) {
            if (isdigit((int)argv[i+1][0]) || isdigit((int)argv[i+1][1])) {
               // first number after keyword is frame at which to freeze
               freeze_at = atoi(argv[++i]);
            }
         }
      } else if (strncmp(argv[i], "-stop", 4) == 0) {
         stop_flow = TRUE;
         if (argc > i+1) {
            if (isdigit((int)argv[i+1][0]) || isdigit((int)argv[i+1][1])) {
               // first number after elevation is frame at which to stop
               stop_at = atoi(argv[++i]);
            }
         }
      } else if (strncmp(argv[i], "-constcn", 8) == 0) {
         courantconst = atof(argv[++i]);

      } else if (strncmp(argv[i], "-tf", 3) == 0) {
         strcpy (tempfilename,argv[++i]);
         use_temp_img = TRUE;
         use_TEMP = TRUE;
      } else if (strncmp(argv[i], "-qf", 3) == 0) {
         strcpy (heatfilename,argv[++i]);
         use_heat_img = TRUE;
         use_TEMP = TRUE;
      } else if (strncmp(argv[i], "-qc", 3) == 0) {
         heat_coeff = atof(argv[++i]);
         use_TEMP = TRUE;
      } else if (strncmp(argv[i], "-tdf", 4) == 0) {
         strcpy (tdfilename,argv[++i]);
         tdlow = atof(argv[++i]);
         tdhigh = atof(argv[++i]);
      } else if (strncmp(argv[i], "-td", 3) == 0) {
         td = atof(argv[++i]);
      } else if (strncmp(argv[i], "-tprint", 3) == 0) {
         print_temp = TRUE;
      } else if (strncmp(argv[i], "-t", 2) == 0) {
         use_TEMP = TRUE;
      } else if (strncmp(argv[i], "-b", 2) == 0) {
         bn = atof(argv[++i]);
         use_strong_strat = FALSE;
      } else if (strncmp(argv[i], "-dr", 3) == 0) {
         dens_ratio = atof(argv[++i]);
         use_strong_strat = TRUE;

      } else if (strncmp(argv[i], "-cf", 3) == 0) {
         strcpy (colorfilename,argv[++i]);
         use_color_img = TRUE;
         use_COLOR = TRUE;
      } else if (strncmp(argv[i], "-cdf", 4) == 0) {
         strcpy (cdfilename,argv[++i]);
         cdlow = atof(argv[++i]);
         cdhigh = atof(argv[++i]);
      } else if (strncmp(argv[i], "-cr", 3) == 0) {
         reread_color = TRUE;
         //fprintf(stderr,"i is %d and argc is %d\n",i,argc);
         //fprintf(stderr,"argv[i] is (%s)\n",argv[i]);
         //fprintf(stderr,"argv[i+1] is (%s)\n",argv[i+1]);
         //fprintf(stderr,"argv[i+1][1] is (%s)\n",argv[i+1][1]);
         //fprintf(stderr,"isdigit argv[i+1][0] is %d\n",isdigit((unsigned char)argv[i+1][0]));
         //fprintf(stderr,"isdigit argv[i+1][1] is %d\n",isdigit((unsigned char)argv[i+1][1]));
         //exit(0);
         if (argc > i+1) {
            if (isdigit(argv[i+1][0]) || isdigit(argv[i+1][1])) {
               overlay_fraction = atof(argv[++i]);
               //if (overlay_fraction < 0.0) {
               //   darkenonly = TRUE;
               //   overlay_fraction *= -1.;
               //}
            }
         }
      } else if (strncmp(argv[i], "-cd", 3) == 0) {
         cd = atof(argv[++i]);
      } else if (strncmp(argv[i], "-c", 2) == 0) {
         use_COLOR = TRUE;

      } else if (strncmp(argv[i], "-df", 3) == 0) {
         strcpy (divfilename,argv[++i]);
         use_div_img = TRUE;
      } else if (strncmp(argv[i], "-div", 3) == 0) {
         use_DF = TRUE;

      } else if (strncmp(argv[i], "-mf", 3) == 0) {
         strcpy (maskfilename,argv[++i]);
         use_mask_img = TRUE;
         use_MASK = TRUE;
      } else if (strncmp(argv[i], "-me", 3) == 0) {
         maskerr = atof(argv[++i]);
      } else if (strncmp(argv[i], "-mprint", 3) == 0) {
         print_mask = TRUE;
      } else if (strncmp(argv[i], "-m", 2) == 0) {
         use_MASK = TRUE;

      } else if (strncmp(argv[i], "-ores", 3) == 0) {
         // NOT IMPLEMENTED
         outscale = atoi(argv[++i]);
      } else if (strncmp(argv[i], "-dt", 3) == 0) {
         dt = atof(argv[++i]);
      } else if (strncmp(argv[i], "-every", 3) == 0) {
         writeevery = atoi(argv[++i]);
         if (writeevery < 1) writeevery = 1;
      } else if (strncmp(argv[i], "-step", 5) == 0) {
         maxstep = atoi(argv[++i]);
      } else if (strncmp(argv[i], "-stam", 5) == 0) {
         isStam = TRUE;
      } else if (strncmp(argv[i], "-fs", 3) == 0) {
         freestream[0] = atof(argv[++i]);
         freestream[1] = atof(argv[++i]);
         xbdry = OPEN;
         ybdry = OPEN;
      } else if (strncmp(argv[i], "-grav", 2) == 0) {
         gravity[0] = atof(argv[++i]);
         gravity[1] = atof(argv[++i]);


      } else if (strncmp(argv[i], "-randvortscale", 6) == 0) {
         randvortscale = atof(argv[++i]);

      } else if (strncmp(argv[i], "-noout", 6) == 0) {
         writeOutput = FALSE;
      } else if (strncmp(argv[i], "-q", 2) == 0) {
         silent = TRUE;
      } else if (strncmp(argv[i], "-8", 2) == 0) {
         use_16bpp = FALSE;
      } else if (strncmp(argv[i], "-16", 3) == 0) {
         use_16bpp = TRUE;
      } else {
         fprintf(stderr,"Unknown option (%s)\n",argv[i]);
         (void) Usage(progname,0);
      }
   }

   // set y max bound (if square, yf=1.0, otherwise scale it as xmax=1.0
   // this is only used in this context in this file, not in the solver!
   yf = (float)(ny-1)/(float)(nx-1);

   // set time step (dt), base it on the maximum courant number (max speed of)
   //  but only if courant number is not set
   if (dt < 0.)
     dt = courant/(nx-1);
   else
     courant = dt*(nx-1);
   if (!silent) fprintf(stdout,"Using dt=%g, or courant number=%g\n",dt,courant);


   // -----------------------------------------------
   // Allocate and initialize

   // velocities first
   u[XV] = allocate_2d_array_f(nx,ny);
   u[YV] = allocate_2d_array_f(nx,ny);

   // one necessary scalar type is vorticity, always at position 0
   a[W2] = allocate_2d_array_f(nx,ny);
   // set coefficient of momentum diffusivity
   sc_diffus[W2] = md;

   // split on solution type
   if (isStam) {

      // both velocities are scalars, as they must be diffused normally
      a[XV] = allocate_2d_array_f(nx,ny);
      // set coefficient of momentum diffusivity
      sc_diffus[XV] = md;

      a[YV] = allocate_2d_array_f(nx,ny);
      sc_diffus[YV] = md;
   }

   // set more scalars; this one: density (or temperature)
   if (use_TEMP) {
      // set thermal diffusivity
      sc_diffus[SF] = td;
      a[SF] = allocate_2d_array_f(nx,ny);
   }

   // this one: divergence
   if (use_DF) {
      // set dilation diffusivity equal to momentum diffusivity,
      // but should really be second viscosity (which is equivalent or greater
      // than momentum diffusivity mu)
      sc_diffus[DF] = md;
      a[DF] = allocate_2d_array_f(nx,ny);
   }

   // this one: color
   if (use_COLOR) {
      // set scalar diffusivity
      sc_diffus[RR] = cd;
      a[RR] = allocate_2d_array_f(nx,ny);
      // now green
      sc_diffus[GG] = cd;
      a[GG] = allocate_2d_array_f(nx,ny);
      // now blue
      sc_diffus[BB] = cd;
      a[BB] = allocate_2d_array_f(nx,ny);

      // save boundary colors
      color_left = allocate_2d_array_f(3,ny);
      color_right = allocate_2d_array_f(3,ny);
      color_top = allocate_2d_array_f(3,nx);
      color_bottom = allocate_2d_array_f(3,nx);
   }

   // variable viscosities

   // variable momentum viscosity
   if (mdlow > 0.0 && mdhigh > 0.0) {
      if (!silent) fprintf(stdout,"Using variable momentum viscosity\n");
      // viscosity diffuses with its own special (constant) diffusivity,
      sc_diffus[MD] = vd;
      a[MD] = allocate_2d_array_f(nx,ny);
      // read grayscale PNG of exactly nx by ny resolution
      read_png(mdfilename,nx,ny,FALSE,FALSE,1.0,FALSE,
         a[MD],mdlow,mdhigh-mdlow,NULL,-1.0,2.0,NULL,-1.0,2.0);
      // and flag the vorticity to use this variable viscosity!
      sc_diffus[W2] = -1.0;
   }

   // allocate space for mask
   if (use_MASK) {
      mask = allocate_2d_array_f(nx,ny);
   } else {
      mask = NULL;
   }

   // allocate space for heat map
   if (use_heat_img) {
      heat = allocate_2d_array_f(nx,ny);
   } else {
      heat = NULL;
   }

   // allocate space for acceleration results
   if (use_strong_strat) {
      acc[0] = allocate_2d_array_f(nx,ny);
      acc[1] = allocate_2d_array_f(nx,ny);
   } else {
      acc[0] = NULL;
      acc[1] = NULL;
   }

   // allocate space for temporary arrays
   for (i=0; i<MAX_SCALARS; i++) {
      if (a[i] != NULL) {
         t[i] = allocate_2d_array_f(nx,ny);
         if (!silent) fprintf(stdout,"Array %d has diffusivity %g\n",i,sc_diffus[i]);
      }
   }

   // -----------------------------------------------
   // Set initial conditions

   // Set vorticity ----------------------------------

   // first, zero it
   for (ix=0; ix<nx; ix++) {
      for (iy=0; iy<ny; iy++) {
         a[W2][ix][iy] = 0.0;
      }
   }

   if (FALSE) {
      // create a wavy shear layer of positive vorticity
      for (ix=0; ix<nx; ix++) {
         px = (float)ix/(float)(nx-1);
         for (iy=0; iy<ny; iy++) {
            py = (float)iy/(float)(nx-1);
            a[W2][ix][iy] = 10.0 - 20.0*fabs((yf/2.0) - py - 0.05*sin(px*2.0*M_PI));
            if (a[W2][ix][iy] < 0.0) a[W2][ix][iy] = 0.0;
         }
      }
   } else if (FALSE) {
      // create three blobs of vorticity (run109/06)
      //add_smooth_circular_blob(nx,ny,xbdry,ybdry,a[W2],0.6,0.2,0.1,10.0);
      //add_smooth_circular_blob(nx,ny,xbdry,ybdry,a[W2],0.4,0.16,0.1,-10.0);
      //add_smooth_circular_blob(nx,ny,xbdry,ybdry,a[W2],0.55,0.5,0.1,-10.0);
      //add_smooth_circular_blob(nx,ny,xbdry,ybdry,a[W2],0.44,0.58,0.1,10.0);
      // run109/07
      add_smooth_circular_blob(nx,ny,xbdry,ybdry,a[W2],0.58,0.37,0.1,10.0);
      add_smooth_circular_blob(nx,ny,xbdry,ybdry,a[W2],0.38,0.33,0.1,-10.0);
      add_smooth_circular_blob(nx,ny,xbdry,ybdry,a[W2],0.53,0.6,0.1,-10.0);
      add_smooth_circular_blob(nx,ny,xbdry,ybdry,a[W2],0.42,0.72,0.1,10.0);
   } else if (FALSE) {
      // create a random field of vorticity
      //scale = 20.;
      scale = 100.;
      for (ix=0; ix<nx; ix++) {
         for (iy=0; iy<ny; iy++) {
            a[W2][ix][iy] = scale*rand()/RAND_MAX - scale*0.5;
         }
      }
   } else if (FALSE) {
      add_smooth_circular_blob(nx,ny,xbdry,ybdry,a[W2],0.6,0.2,0.1,10.0);
      add_smooth_circular_blob(nx,ny,xbdry,ybdry,a[W2],0.4,0.16,0.1,-10.0);
   } else if (FALSE) {
      add_smooth_circular_blob(nx,ny,xbdry,ybdry,a[W2],0.57,0.37,0.05,50.0);
      add_smooth_circular_blob(nx,ny,xbdry,ybdry,a[W2],0.41,0.34,0.05,-50.0);
   } else if (FALSE) {
      // make the vorticity defined in Minion and Brown, JCP 138
      for (ix=0; ix<nx; ix++) {
         px = (float)ix/(float)(nx-1);
         for (iy=0; iy<ny; iy++) {
            py = (float)iy/(float)(nx-1);
            if (py > 0.5) {
               a[W2][ix][iy] = 2.*M_PI*0.05*cos(2.*M_PI*(px+0.25))
                             + 80./pow(0.5*(exp(80.*(0.75-py))+exp(80.*(py-0.75))),2);
            } else {
               a[W2][ix][iy] = 2.*M_PI*0.05*cos(2.*M_PI*(px+0.25))
                             - 80./pow(0.5*(exp(80.*(0.25-py))+exp(80.*(py-0.25))),2);
            }
         }
      }
   } else if (FALSE) {
      // make the vorticity defined in Koumoutsakos, JCP, 1997, type I
      for (ix=0; ix<nx; ix++) {
         px = (float)ix/(float)(nx-1);
         for (iy=0; iy<ny; iy++) {
            py = (float)iy/(float)(nx-1);
            temp = sqrt(pow(px-0.5,2)+0.25*pow(py-0.5,2))/0.1;
            temp2 = 20.*(1.-exp(-exp(1./(temp-1.))*2.56085/temp));
            if (temp2 < 0.) temp2 = 0.;
            if (temp2 > 20.) temp2 = 20.;
            if (temp < 1.) a[W2][ix][iy] = temp2;
            //if (temp < 1.) fprintf(stderr,"%d %d %g\n",ix,iy,temp2);
            //if (isnan(temp2)) fprintf(stderr,"%d %d %g\n",ix,iy,temp2);
         }
      }
   } else if (FALSE) {
      // create a sinusoidal interface
      //thickness = 0.001;
      thickness = 2.0/(float)nx;
      for (ix=0; ix<nx; ix++) {
         px = (ix - (float)nx/2.)/(float)nx;
         yf = 0.02*sin(px*9.*M_PI);
         for (iy=0; iy<ny; iy++) {
            py = (iy - (float)ny/2.)/(float)ny;
            // now, px,py is [-0.5,0.5]
            if (fabs(yf-py) < thickness) {
               // point is very close to the curve
               temp = 0.5+0.5*cos(M_PI*(yf-py)/thickness);
               //fprintf(stderr,"%g %g\n",fabs(yf-py),temp);
               temp += 0.02*rand()/RAND_MAX - 0.01;
               temp *= 200.*cos(px*9.*M_PI);
               a[W2][ix][iy] = temp;
            //} else {
               //temp = 0.;
            }
            //fprintf(stderr,"%d %d  %g %g  %g  %g\n",ix,iy,px,py,yf,temp);
         }
      }
   }
   if (use_vort_img) {
      // read grayscale PNG of exactly nx by ny resolution
      // the funny min bounds are to allow value of 127 to become 0.0
      read_png(vortfilename,nx,ny,FALSE,
               FALSE,1.0,FALSE,
               a[W2],-254.*vortscale/255.,2.0*vortscale,
               NULL, -254.*vortscale/255.,2.0*vortscale,
               NULL, -254.*vortscale/255.,2.0*vortscale);
   }
   if (randvortscale > 0.0) {
      // create a random field of vorticity
      //scale = 20.;
      scale = randvortscale;
      for (ix=0; ix<nx; ix++) {
         for (iy=0; iy<ny; iy++) {
            a[W2][ix][iy] += scale*(rand()/(float)RAND_MAX - 0.5);
         }
      }
   }

   // Set wall velocity BCs --------------------------
   // do not use this, it doesn't work
   //int iwall = -1;
   //funcbndyc (&iwall, (float*)NULL, (float*)NULL, &wallvel[0]);

   // Set flow mask ----------------------------------

   if (use_MASK) {
      // define a mask manually
      // blank the mask
      for (ix=0; ix<nx; ix++) {
         for (iy=0; iy<ny; iy++) {
            // use 0.0 for all solid
            mask[ix][iy] = 1.0;
         }
      }
      if (FALSE) {
         // create a circular mask at cx,cy
         cx = 0.3;
         cy = 0.5;
         rad = 0.15;
         //temp2 = 0.0025;
         temp2 = 2.0/(float)nx;
         rad -= temp2/2.;
         for (ix=0; ix<nx; ix++) {
            px = (float)ix/(float)(nx-1);
            for (iy=0; iy<ny; iy++) {
               py = (float)iy/(float)(nx-1);
               temp = sqrt(pow(px-cx,2)+pow(py-cy,2));
               if (temp < rad) {
                  // use 1.0 for a fluid hole
                  mask[ix][iy] = 0.0;
               } else if (temp < rad+temp2) {
                  // use this first construct for a central hole
                  //mask[ix][iy] = cos(10.*(temp-0.1)*M_PI+0.5);
                  // use this construct for a central solid
                  mask[ix][iy] = 0.5-0.5*cos((temp-rad)*M_PI/temp2);
               }
            }
         }
      }
      if (FALSE) {
         // create a solar tower
         cx = 0.5;	// x center
         cy = 0.0;	// y center
         minval = 0.01; // inner surface
         maxval = 0.06; // outer surface
         temp2 = 0.002;	// thickness
         scale = 0.3;	// length
         for (ix=0; ix<nx; ix++) {
           px = fabs((float)ix/(float)(nx-1) - cx);
           if (px < scale) {
            for (iy=0; iy<ny; iy++) {
              py = (float)iy/(float)(nx-1) - cy;
              if (py > 0. && py < scale) {
              //fprintf(stderr,"%d %d %g %g %g\n",ix,iy,px,py,px*py);
               temp = px*py;
               if (temp > minval && temp < maxval) {
                  // use 0.0 for a solid
                  mask[ix][iy] = 0.0;
               } else if (temp > maxval && temp < maxval+temp2) {
                  // use this construct for a central solid
                  mask[ix][iy] = (temp-maxval)/temp2;
               } else if (temp < minval && temp > minval-temp2) {
                  // use this construct for a central solid
                  mask[ix][iy] = (minval-temp)/temp2;
               }
              }
            }
           }
         }
      }
      // read in the mask file
      if (use_mask_img) {
         // read grayscale PNG of exactly nx by ny resolution
         read_png(maskfilename,nx,ny,FALSE,FALSE,1.0,FALSE,
            mask,0.0,1.0,NULL,0.0,1.0,NULL,0.0,1.0);
      }
      // normalize mask, and set 1.0=solid, 0.0=open (input is opposite)
      maxval = -9.9e+9;
      minval = 9.9e+9;
      for (ix=0; ix<nx; ix++) {
         for (iy=0; iy<ny; iy++) {
            if (mask[ix][iy] > maxval) maxval = mask[ix][iy];
            if (mask[ix][iy] < minval) minval = mask[ix][iy];
         }
      }
      temp = maxval-minval;
      if (temp > 1.e-5) {
      for (ix=0; ix<nx; ix++) {
         for (iy=0; iy<ny; iy++) {
            mask[ix][iy] = (maxval-mask[ix][iy])/temp;
            if (mask[ix][iy] < 0.0) mask[ix][iy] = 0.0;
            if (mask[ix][iy] > 1.0) mask[ix][iy] = 1.0;
         }
      }
      }
   }

   // optionally generate repeatedly-overlaid mask
   //(void) overlay_mask (nx, ny, mask);

   // Initial velocity solve -------------------------

   // if you want to show velocity on the first step, solve for it here
   find_vels_2d (silent,step,isStam,nx,ny,xbdry,ybdry,freestream,wallvel,u[XV],u[YV],a[W2],use_MASK,mask,maskerr);

   // find vmax (might need this)
   vmax = 0.;
   for (i=0; i<nx; i++) {
      for (j=0; j<ny; j++) {
         velsq = u[XV][i][j]*u[XV][i][j] + u[YV][i][j]*u[YV][i][j];
         if (velsq > vmax) vmax = velsq;
      }
   }
   vmax = sqrt(vmax);
   // initialize the vmax averaging array
   vmax = courantconst / (dt*(nx+1));
   for (i=0; i<VMAXAVG; i++) ccnvmax[i] = vmax;


   // set color ---------------------------------------

   // load in an image containing colors to grab
   if (FALSE) {

      // the color image from which to grab colors
      sprintf(colorsrcfilename,"color_source.png");
      cnx = 300;
      cny = 169;

      // allocate space and read image
      c[0] = allocate_2d_array_f(cnx,cny);
      c[1] = allocate_2d_array_f(cnx,cny);
      c[2] = allocate_2d_array_f(cnx,cny);
      read_png (colorsrcfilename,cnx,cny,TRUE,FALSE,1.0,FALSE,
               c[0],0.0,1.0,c[1],0.0,1.0,c[2],0.0,1.0);

      // later on, grab colors with
      //(void) get_random_color (c,cnx,cny,thisc);
   }

   // initialize an array of blocks to add during the run
   if (FALSE) {
      populate_block_array(nx,ny);
   }

   // first, find what index the scalar uses
   if (use_COLOR) {
      // zero the arrays
      //for (ix=0; ix<nx; ix++) for (iy=0; iy<ny; iy++) a[RR][ix][iy] = 0.0;
      //for (ix=0; ix<nx; ix++) for (iy=0; iy<ny; iy++) a[GG][ix][iy] = 0.0;
      //for (ix=0; ix<nx; ix++) for (iy=0; iy<ny; iy++) a[BB][ix][iy] = 0.0;
      //for (ix=0; ix<nx; ix++) for (iy=0; iy<ny; iy++) a[RR][ix][iy] = 1.0;
      //for (ix=0; ix<nx; ix++) for (iy=0; iy<ny; iy++) a[GG][ix][iy] = 1.0;
      //for (ix=0; ix<nx; ix++) for (iy=0; iy<ny; iy++) a[BB][ix][iy] = 1.0;
      for (ix=0; ix<nx; ix++) for (iy=0; iy<ny; iy++) a[RR][ix][iy] = 0.10546875;
      for (ix=0; ix<nx; ix++) for (iy=0; iy<ny; iy++) a[GG][ix][iy] = 0.15234375;
      for (ix=0; ix<nx; ix++) for (iy=0; iy<ny; iy++) a[BB][ix][iy] = 0.05859375;
      if (use_color_img) {
         // read color PNG of exactly nx by ny resolution into
         //   the given fields
         read_png (colorfilename,nx,ny,TRUE,FALSE,1.0,FALSE,
            a[RR],0.0,1.0,a[GG],0.0,1.0,a[BB],0.0,1.0);
         //read_png (colorfilename,nx,ny,TRUE,FALSE,1.0,FALSE,
         //   a[RR],-3.0,4.0,a[GG],-3.0,4.0,a[BB],-3.0,4.0);
      }
      if (FALSE) {
         // left is red, right is green
         for (ix=0; ix<nx/2; ix++) {
            for (iy=0; iy<ny; iy++) {
               a[RR][ix][iy] = 1.0;
               a[GG][ix][iy] = 0.91372549;
               a[BB][ix][iy] = 0.168627451;
            }
         }
         for (ix=nx/2; ix<nx; ix++) {
            for (iy=0; iy<ny; iy++) {
               a[RR][ix][iy] = 0.1171875;
               a[GG][ix][iy] = 0.0546875;
               a[BB][ix][iy] = 0.7265625;
            }
         }
      }
      if (FALSE) {
         // first, force background to one color
         for (ix=0; ix<nx; ix++) for (iy=0; iy<ny; iy++) a[RR][ix][iy] = 235./256.;
         for (ix=0; ix<nx; ix++) for (iy=0; iy<ny; iy++) a[GG][ix][iy] = 245./256.;
         for (ix=0; ix<nx; ix++) for (iy=0; iy<ny; iy++) a[BB][ix][iy] = 253./256.;
         // make splotches of colored spots
         // startx,starty, endx,endy, radius,cutoff, r,g,b, numspots, nx,ny,ra,ga,ba
         // blue splat
         //(void) paint_splat (0.27,0.3, 0.67,0.5, 0.1,2.e+4,
         //                    0.1171875,0.0546875,0.7265625, 1000,
         //                    nx,ny,a[RR],a[GG],a[BB]);
         // variety of splats
         (void) paint_splat (-1.5,0.0, 1.3,1.0, 0.5,0.9e+3,
                             141./256.,206./256.,228./256., 1000,
                             nx,ny,a[RR],a[GG],a[BB]);
         (void) paint_splat (-0.5,-0.5, 0.5,2.0, 0.4,1.0e+3,
                             191./256.,226./256.,248./256., 1000,
                             nx,ny,a[RR],a[GG],a[BB]);
         (void) paint_splat (2.1,-0.3, 0.9,0.25, 0.3,1.7e+3,
                             141./256.,206./256.,228./256., 1000,
                             nx,ny,a[RR],a[GG],a[BB]);
         (void) paint_splat (-0.2,-0.4, 0.49,0.5, 0.25,2.e+3,
                             124./256.,155./256.,193./256., 1000,
                             nx,ny,a[RR],a[GG],a[BB]);
      }
      if (FALSE) {
         // scale input vorticity by pixel brightness
         for (ix=0; ix<nx; ix++) {
            for (iy=0; iy<ny; iy++) {
               float brite = 0.3*a[RR][ix][iy] + 0.6*a[GG][ix][iy] + 0.1*a[BB][ix][iy];
               a[W2][ix][iy] *= 1.0-brite;
            }
         }
      }

      // save the edge colors
      for (ix=0; ix<nx; ix++) {
         color_bottom[0][ix] = a[RR][ix][0];
         color_bottom[1][ix] = a[GG][ix][0];
         color_bottom[2][ix] = a[BB][ix][0];
         color_top[0][ix] = a[RR][ix][ny-1];
         color_top[1][ix] = a[GG][ix][ny-1];
         color_top[2][ix] = a[BB][ix][ny-1];
      }
      for (iy=0; iy<ny; iy++) {
         color_left[0][iy] = a[RR][0][iy];
         color_left[1][iy] = a[GG][0][iy];
         color_left[2][iy] = a[BB][0][iy];
         color_right[0][iy] = a[RR][nx-1][iy];
         color_right[1][iy] = a[GG][nx-1][iy];
         color_right[2][iy] = a[BB][nx-1][iy];
      }
   }

   // Set scalar/temperature ----------------------------------

   // first, find what index the scalar uses
   if (use_TEMP) {

      // zero the array
      for (ix=0; ix<nx; ix++) {
         for (iy=0; iy<ny; iy++) {
            a[SF][ix][iy] = 0.0;
         }
      }

      if (use_temp_img) {
         // read grayscale PNG of exactly nx by ny resolution
         //read_png(tempfilename,nx,ny,FALSE,TRUE,1.0,FALSE,
         read_png(tempfilename,nx,ny,FALSE,FALSE,1.0,FALSE,
            a[SF],-1.0,2.0,NULL,0.0,1.0,NULL,0.0,1.0);
      }

      if (use_heat_img) {
         // read grayscale PNG of exactly nx by ny resolution
         read_png(heatfilename,nx,ny,FALSE,FALSE,1.0,FALSE,
            heat,-1.0,2.0,NULL,0.0,1.0,NULL,0.0,1.0);

         for (ix=0; ix<nx; ix++) {
            for (iy=0; iy<ny; iy++) {
               a[SF][ix][iy] += heat[ix][iy];
               if (a[SF][ix][iy] > 1.0) a[SF][ix][iy] = 1.0;
               if (a[SF][ix][iy] < -1.0) a[SF][ix][iy] = -1.0;
            }
         }
      }

      if (FALSE) {
         // create a sharp blob of scalar
         for (ix=0; ix<nx; ix++) {
            px = (float)ix/(float)(nx-1);
            for (iy=0; iy<ny; iy++) {
               py = (float)iy/(float)(nx-1);
               if (sqrt(pow(px-0.4,2)+pow(py-0.15,2)) < 0.1)
                  a[SF][ix][iy] = 10.0;
            }
         }
      }
      if (FALSE) {
         // create a smooth blob of scalar
         // add_smooth_circular_blob(nx,ny,xbdry,ybdry,a[SF],0.7,0.15,0.1,1.0);
         add_smooth_circular_blob(nx,ny,xbdry,ybdry,a[SF],0.4,0.0,0.07,-1.0);
      }
      if (FALSE) {
         // create bands of scalar
         for (ix=0; ix<nx; ix++) {
            for (iy=0; iy<ny; iy++) {
               //a[SF][ix][iy] = sin((double)(ix*0.2)) +
               //                 0.2*rand()/RAND_MAX - 0.1;
               //a[SF][ix][iy] = 5.*sin((double)(ix*0.2));
               // run109/06
               //a[SF][ix][iy] = 5.*sin(410.*(double)(ix)/nx);
               // run109/07
               a[SF][ix][iy] = 0.5 + 0.5*sin(180.*(double)(ix)/nx);
            }
         }
      }
      if (FALSE) {
         // create square bands of scalar
         for (ix=0; ix<nx; ix++) {
            for (iy=0; iy<ny; iy++) {
               //a[SF][ix][iy] = sin((double)(ix*0.2)) +
               //                 0.2*rand()/RAND_MAX - 0.1;
               px = fabs(ix - (float)nx/2.);
               py = fabs(iy - (float)ny/2.);
               if (px > py) temp = px;
               else temp = py;
               if (temp > nx/3.1) temp = -M_PI/2.;
               else temp *= 0.3;
               a[SF][ix][iy] = 5.*sin(temp);
            }
         }
      }
      if (FALSE) {
         // create a RTI with a sinusoidal interface
         //thickness = 0.001;
         thickness = 2.0/(float)nx;
         //thickness = 1.0/(float)nx;
         for (ix=0; ix<nx; ix++) {
            px = (ix - (float)nx/2.)/(float)nx;
            //yf = 0.02*sin(px*9.*M_PI);
            yf = 0.02*sin(px*8.*M_PI);
            for (iy=0; iy<ny; iy++) {
               //a[SF][ix][iy] = sin((double)(ix*0.2)) +
               //                 0.2*rand()/RAND_MAX - 0.1;
               py = (iy - (float)ny/2.)/(float)ny;
               // now, px,py is [-0.5,0.5]
               if (yf-py > thickness) {
                  // point is below the curve
                  temp = 0.;
               } else if (yf-py < -thickness) {
                  // point is below the curve
                  temp = 1.;
               } else {
                  // point is very close to the curve
                  temp = 0.5-0.5*sin(0.5*M_PI*(yf-py)/thickness);
                  temp += 0.02*rand()/RAND_MAX - 0.01;
               }
               //fprintf(stderr,"%d %d  %g %g  %g  %g\n",ix,iy,px,py,yf,temp);
               a[SF][ix][iy] = temp;
            }
         }
      }
      if (FALSE) {
         // create another sinusoidal interface (bottom one)
         //thickness = 0.001;
         thickness = 2.0/(float)nx;
         for (ix=0; ix<nx; ix++) {
            px = (ix - (float)nx/2.)/(float)nx;
            yf = 0.02*sin(px*9.*M_PI) - 0.02;
            for (iy=0; iy<ny; iy++) {
               //a[SF][ix][iy] = sin((double)(ix*0.2)) +
               //                 0.2*rand()/RAND_MAX - 0.1;
               py = (iy - (float)ny/2.)/(float)ny;
               // now, px,py is [-0.5,0.5]
               if (yf-py > thickness) {
                  // point is below the curve
                  temp = 1.;
               } else if (yf-py < -thickness) {
                  // point is above the curve
                  temp = 0.;
               } else {
                  // point is very close to the curve
                  temp = 0.5+0.5*sin(0.5*M_PI*(yf-py)/thickness);
                  temp += 0.02*rand()/RAND_MAX - 0.01;
               }
               //fprintf(stderr,"%d %d  %g %g  %g  %g\n",ix,iy,px,py,yf,temp);
               a[SF][ix][iy] += temp;
            }
         }
      }
      if (FALSE) {
        // create a random field of scalar
        scale = 10.;
        for (ix=0; ix<nx; ix++) {
          for (iy=0; iy<ny; iy++) {
            a[SF][ix][iy] += scale*rand()/RAND_MAX - 0.5*scale + 0.5;
          }
        }
      }
      if (FALSE && use_COLOR) {
         // use the brightness component to indicate density!
         for (ix=0; ix<nx; ix++) {
            for (iy=0; iy<ny; iy++) {
               a[SF][ix][iy] = 1.0 - 0.6*a[GG][ix][iy] -
                                0.3*a[RR][ix][iy] -
                                0.1*a[BB][ix][iy];
            }
         }
      }
   }

   // -----------------------------------------------
   // Iterate through time

   step = 0;
   simtime = 0.0;
   totcputime = 0.0;
   while (TRUE) {

      // fprintf(stdout,"\nBegin step %d\n",step);

      // ----------------------------
      // write output
      if (writeOutput && step%writeevery == 0) {
      //#pragma omp parallel sections private(outfileroot)
      {
         //#pragma omp section
         if (print_vort) {
            sprintf(outfileroot,"vort_%06d",step);
            //write_png (outfileroot,nx,ny,FALSE,FALSE,a[W2],-1.0,2.0,
            //write_png (outfileroot,nx,ny,FALSE,FALSE,a[W2],-10.0,20.0,
            //write_png (outfileroot,nx,ny,FALSE,FALSE,a[W2],-250.0,500.0,
            //write_png (outfileroot,nx,ny,FALSE,FALSE,a[W2],-100.0,200.0,
            //write_png (outfileroot,nx,ny,FALSE,FALSE,a[W2],-10.0,20.0,
            write_png (outfileroot,nx,ny,FALSE,use_16bpp,
                       a[W2],-vortscale,2.*vortscale,
                       NULL,0.0,1.0,
                       NULL,0.0,1.0);
         }
         //#pragma omp section
         if (print_mask && use_MASK) {
            sprintf(outfileroot,"mask_%06d",step);
            write_png (outfileroot,nx,ny,FALSE,FALSE,
                       mask,0.0,1.0,
                       NULL,0.0,1.0,
                       NULL,0.0,1.0);
         }
         //#pragma omp section
         if (print_temp) {
            sprintf(outfileroot,"temp_%06d",step);
            write_png (outfileroot,nx,ny,FALSE,FALSE,
                       a[SF],-1.0,2.0,
                       NULL,0.0,1.0,
                       NULL,0.0,1.0);
         }
         //#pragma omp section
         if (use_COLOR) {
            sprintf(outfileroot,"out_%06d",step);
            write_png (outfileroot,nx,ny,TRUE,use_16bpp,
                       a[RR],0.0,1.0,
                       a[GG],0.0,1.0,
                       a[BB],0.0,1.0);
         }
         //#pragma omp section
         if (print_vel) {
            sprintf(outfileroot,"vel_%06d",step);
            write_png (outfileroot,nx,ny,TRUE,use_16bpp,
                       u[XV],-velscale,2.*velscale,
                       u[YV],-velscale,2.*velscale,
                       a[W2],-vortscale,2.*vortscale);
         }
         //#pragma omp section
         if (print_mu) {
            sprintf(outfileroot,"mu_%06d",step);
            write_png (outfileroot,nx,ny,FALSE,use_16bpp,
                       a[MD],mdlow,mdhigh-mdlow,
                       NULL,0.0,1.0,
                       NULL,0.0,1.0);
         }
      } // OMP sections
      } // if writeOutput && step%writeevery == 0

      // check end conditions
      if (step >= maxstep) break;

      // accept input (only if interactive scheme)

      // adjust diffusion coefficients - decay past step 300
      if (FALSE && step > 300) {
        for (i=0; i<MAX_SCALARS; i++) {
          // this will work for the negative diffusivity trigger, too
          sc_diffus[i] *= 1.1;
        }
      }

      // if we want a constant cn time step, compute it here
      if (courantconst > 0.0) {
         // average vmax over a number of steps
         ccnidx = (ccnidx+1)%VMAXAVG;
         ccnvmax[ccnidx] = vmax;
         // if any were zero, reset them
         for (i=0; i<VMAXAVG; i++) if (ccnvmax[i] < 1.e-6) ccnvmax[i] = vmax;
         // find new thing
         vmax = 0.0;
         for (i=0; i<VMAXAVG; i++) vmax += ccnvmax[i];
         vmax /= (float)VMAXAVG;
         // recalculate dt
         float newdt = (courantconst/vmax)/(nx+1);
         // do not let dt change too much!
         if (newdt > 2.0*dt) dt *= 2.0;
         else if (newdt < 0.5*dt) dt *= 0.5;
         else dt = newdt;
         if (!silent) fprintf(stderr,"  changing dt to %g\n",dt);
      }

      // set the timer
      last_tics = clock();

      // should we test for freeze?
      if (freeze_flow) {
         if (step > freeze_at) {
            recalc_vel = FALSE;
            move_colors = TRUE;
         } else {
            recalc_vel = TRUE;
            move_colors = FALSE;
         }
      }

      // should we test for stop?
      if (stop_flow) {
         if (step > stop_at) {
            recalc_vel = FALSE;
            move_colors = FALSE;
         } else {
            recalc_vel = TRUE;
            move_colors = TRUE;
         }
      }

      // take one computational step forward in time
      int numsubsteps = 1;
      if (use_MASK) numsubsteps = 1 + (int)(vmax*(nx+1)*dt/10.);
      for (int istep=0; istep<numsubsteps; istep++) {
         effective_re = step_forward_2d (silent,step,isStam,4,
                           nx,ny,xbdry,ybdry,freestream,wallvel,
                           recalc_vel,move_colors,
                           u,a,t,use_MASK,mask,maskerr,sc_diffus,
                           gravtype,gravity,(dt/(float)numsubsteps),
                           use_strong_strat,bn,dens_ratio,acc);
      }

      // read in the image again and overlay it!
      if (use_COLOR && reread_color) {
         // read color PNG of exactly nx by ny resolution into
         //   the given fields and overlay it!
         if (use_color_img && reread_color) {
            //read_png(colorfilename,nx,ny,TRUE,TRUE,0.01,TRUE,
            //   a[RR],0.0,1.0,a[GG],0.0,1.0,a[BB],0.0,1.0);
            read_png(colorfilename,nx,ny,TRUE,
                     TRUE,overlay_fraction,darkenonly,
                     a[RR],0.0,1.0,a[GG],0.0,1.0,a[BB],0.0,1.0);
            // use the brightness component to indicate density!
            //for (ix=0; ix<nx; ix++) {
            //   for (iy=0; iy<ny; iy++) {
            //      a[SF][ix][iy] = 1.0 - 0.6*a[GG][ix][iy] -
            //                       0.3*a[RR][ix][iy] -
            //                       0.1*a[BB][ix][iy];
            //   }
            //}
         }
      }

      // read in the vorticity image and overlay it
      if (recalc_vel && reread_vort && use_vort_img) {
         // read grayscale PNG of exactly nx by ny resolution
         // the 2 is to force the new data to be simply added to the old
         // the funny min bounds are to allow value of 127 to become 0.0
         read_png(vortfilename,nx,ny,FALSE,
                  2,1.0,FALSE,
                  a[W2],-254.*vortscale/255.,2.0*vortscale,
                  NULL, -254.*vortscale/255.,2.0*vortscale,
                  NULL, -254.*vortscale/255.,2.0*vortscale);
      }

      if (FALSE) {
         // open or close blocks in the mask
         (void) update_mask_with_blocks_2 (mask, a, nx, ny, simtime, dt);
         // optionally generate repeatedly-overlaid mask
         (void) overlay_mask (nx, ny, mask);
      }

      if (FALSE && step%50 == 0 && step < 1001) {
         printf("Adding splotch %d\n",step);
         // grab a color and splat the paint
         (void) get_random_color (c,cnx,cny,thisc);
         (void) paint_splat (-1.0+3.0*rand()/(float)RAND_MAX, -1.0+3.0*rand()/(float)RAND_MAX,
                             -1.0+3.0*rand()/(float)RAND_MAX, -1.0+3.0*rand()/(float)RAND_MAX,
                             0.1+0.5*rand()/(float)RAND_MAX,
                             (2.e+3)*(0.5+rand()/(float)RAND_MAX),
                             thisc[0],thisc[1],thisc[2], 1000,
                             nx,ny,a[RR],a[GG],a[BB]);
         // reset scalar temperature
         for (ix=0; ix<nx; ix++) {
            for (iy=0; iy<ny; iy++) {
               a[SF][ix][iy] = 1.0 - 0.6*a[GG][ix][iy] -
                                0.3*a[RR][ix][iy] -
                                0.1*a[BB][ix][iy];
            }
         }
      }

      if (FALSE && use_MASK) {
         // allocate space
         if (first_time) {
            shear = allocate_2d_array_f(nx,ny);
         }

         // adjust color value based on shear present in flow
         (void) find_shear_magnitude (nx, ny, xbdry, ybdry, u[XV],1.0, u[YV],1.0, shear);

         float minshear = 9.9e+9;
         float maxshear = -9.9e+9;
         float meanshear = 0.0;

         for (ix=0; ix<nx; ix++) {
            for (iy=0; iy<ny; iy++) {
               if (shear[ix][iy] < minshear) minshear = shear[ix][iy];
               if (shear[ix][iy] > maxshear) maxshear = shear[ix][iy];

               // adjust color
               //float shearscale = 1.0 - 0.002*(fabs(shear[ix][iy]) - 1.0);
               //a[RR][ix][iy] *= shearscale;
               //a[GG][ix][iy] *= shearscale;
               //a[BB][ix][iy] *= shearscale;

               // or adjust mask
               mask[ix][iy] += 1.e-5 * fabs(shear[ix][iy]);
               if (mask[ix][iy] > 1.0) mask[ix][iy] = 1.0;
               if (mask[ix][iy] < 0.0) mask[ix][iy] = 0.0;
            }
         }
         printf("  shear range: %g to %g, mean abs is %g\n", minshear, maxshear, meanshear/(float)(nx*ny));

         if (TRUE) {
            sprintf(outfileroot,"shear_%06d",step);
            write_png (outfileroot,nx,ny,FALSE,FALSE,
                       shear,-100.0,200.0,
                       NULL,0.0,1.0,
                       NULL,0.0,1.0);
         }
      }

      // use vorticity and mask to modify mask
      if (FALSE && use_MASK) {
         // if vorticity is less than this, deposit "mask"
         const float depos_thresh = 100.0;
         // if vorticity is greater than this, erode "mask"
         const float erode_thresh = 200.0;

         float maxvort, thisvort, minmask, maxmask, thismask;
         int ixx, iyy;
         for (ix=1; ix<nx-1; ix++) {
            for (iy=1; iy<ny-1; iy++) {
               // find min/max of neighbors
               minmask = 1.0;
               maxmask = 0.0;
               for (ixx=ix-1; ixx<ix+2; ixx++) {
                  for (iyy=iy-1; iyy<iy+2; iyy++) {
                     thismask = mask[ixx][iyy];
                     if (thismask < minmask) minmask = thismask;
                     if (thismask > maxmask) maxmask = thismask;
                  }
               }
               maxvort = 0.0;
               for (ixx=ix-1; ixx<ix+2; ixx++) {
                  for (iyy=iy-1; iyy<iy+2; iyy++) {
                     thisvort = fabs(a[W2][ixx][iyy]);
                     if (thisvort > maxvort) maxvort = thisvort;
                  }
               }
               // note: mask = 0 means that there is solid stuff there, 1.0 means that it is wide open
               // now test vs. thresholds
               if (maxvort < depos_thresh) {
                  // to deposit, there needs to be something nearby
                  if (minmask < 0.1) {
                     // deposit mask
                     mask[ix][iy] -= 1.e-3 * (depos_thresh-maxvort);
                     // and the new mask value can't be lower than the nearby minimum
                     if (mask[ix][iy] < minmask) mask[ix][iy] = minmask;
                  }
               }
               if (maxvort > erode_thresh) {
                  // to erode, there needs to be some fluid nearby (very likely)
                  if (maxmask > 0.9) {
                     // erode mask
                     mask[ix][iy] += 1.e-3 * (maxvort - erode_thresh);
                     // and the new mask value can't be higher than the nearby maximum
                     if (mask[ix][iy] > maxmask) mask[ix][iy] = maxmask;
                  }
               }

               // maintain mask bounds
               if (mask[ix][iy] > 1.0) mask[ix][iy] = 1.0;
               if (mask[ix][iy] < 0.0) mask[ix][iy] = 0.0;
            }
         }
      }

      // use velocity and mask to modify mask
      if (FALSE && use_MASK) {
         // if velocity magnitude is less than this, deposit "mask"
         const float depos_thresh = 0.2;
         // if velocity magnitude is greater than this, erode "mask"
         const float erode_thresh = 0.2;

         float maxvel, thisvel, minmask, maxmask, thismask;
         int ixx, iyy;

         for (ix=1; ix<nx-1; ix++) {
            for (iy=1; iy<ny-1; iy++) {
               // find min/max of neighbors
               minmask = 1.0;
               maxmask = 0.0;
               for (ixx=ix-1; ixx<ix+2; ixx++) {
                  for (iyy=iy-1; iyy<iy+2; iyy++) {
                     thismask = mask[ixx][iyy];
                     if (thismask < minmask) minmask = thismask;
                     if (thismask > maxmask) maxmask = thismask;
                  }
               }
               maxvel = 0.0;
               for (ixx=ix-1; ixx<ix+2; ixx++) {
                  for (iyy=iy-1; iyy<iy+2; iyy++) {
                     thisvel = sqrt(pow(u[XV][ixx][iyy],2) + pow(u[YV][ixx][iyy],2));
                     if (thisvel > maxvel) maxvel = thisvel;
                  }
               }
               // note: mask = 0 means that there is solid stuff there, 1.0 means that it is wide open
               // now test vs. thresholds
               if (maxvel < depos_thresh) {
                  // to deposit, there needs to be something nearby
                  //if (minmask < 0.1) {
                     // deposit mask
                     mask[ix][iy] -= 10.0*dt * (depos_thresh - maxvel);
                     // and the new mask value can't be lower than the nearby minimum
                     if (mask[ix][iy] < minmask) mask[ix][iy] = minmask;
                  //}
               }
               if (maxvel > erode_thresh) {
                  // to erode, there needs to be some fluid nearby (very likely)
                  //if (maxmask > 0.9) {
                     // erode mask
                     mask[ix][iy] += 10.0*dt * (maxvel - erode_thresh);
                     // and the new mask value can't be higher than the nearby maximum
                     if (mask[ix][iy] > maxmask) mask[ix][iy] = maxmask;
                  //}
               }

               // maintain mask bounds
               if (mask[ix][iy] > 1.0) mask[ix][iy] = 1.0;
               if (mask[ix][iy] < 0.0) mask[ix][iy] = 0.0;
            }
         }
      }

      // overlay the heat map
      if (use_heat_img) {
         for (ix=0; ix<nx; ix++) {
            for (iy=0; iy<ny; iy++) {
               a[SF][ix][iy] += dt * heat_coeff * heat[ix][iy];
               if (a[SF][ix][iy] > 1.0) a[SF][ix][iy] = 1.0;
               if (a[SF][ix][iy] < -1.0) a[SF][ix][iy] = -1.0;
            }
         }
      }

      // overlay the edge color input
      if (use_COLOR && overlay_color_tb) {
         int nthick = 6;
         for (iy=0; iy<nthick; iy++) {
            float tfac = 1.0 / (float)(iy+1);
            for (ix=0; ix<nx; ix++) {
               a[RR][ix][iy] = tfac*color_bottom[0][ix] + (1.0-tfac)*a[RR][ix][iy];
               a[GG][ix][iy] = tfac*color_bottom[1][ix] + (1.0-tfac)*a[GG][ix][iy];
               a[BB][ix][iy] = tfac*color_bottom[2][ix] + (1.0-tfac)*a[BB][ix][iy];
               a[RR][ix][ny-1-iy] = tfac*color_top[0][ix] + (1.0-tfac)*a[RR][ix][ny-1-iy];
               a[GG][ix][ny-1-iy] = tfac*color_top[1][ix] + (1.0-tfac)*a[GG][ix][ny-1-iy];
               a[BB][ix][ny-1-iy] = tfac*color_top[2][ix] + (1.0-tfac)*a[BB][ix][ny-1-iy];
            }
         }
      }
      if (use_COLOR && overlay_color_lr) {
         int nthick = 6;
         for (ix=0; ix<nthick; ix++) {
            float tfac = 1.0 / (float)(ix+1);
            for (iy=0; iy<ny; iy++) {
               a[RR][ix][iy] = tfac*color_left[0][iy] + (1.0-tfac)*a[RR][ix][iy];
               a[GG][ix][iy] = tfac*color_left[1][iy] + (1.0-tfac)*a[GG][ix][iy];
               a[BB][ix][iy] = tfac*color_left[2][iy] + (1.0-tfac)*a[BB][ix][iy];
               a[RR][nx-1-ix][iy] = tfac*color_right[0][iy] + (1.0-tfac)*a[RR][nx-1-ix][iy];
               a[GG][nx-1-ix][iy] = tfac*color_right[1][iy] + (1.0-tfac)*a[GG][nx-1-ix][iy];
               a[BB][nx-1-ix][iy] = tfac*color_right[2][iy] + (1.0-tfac)*a[BB][nx-1-ix][iy];
            }
         }
      }

      // calculate the total time elapsed, and the time for this last step
      tics = clock();
      cputime = ((float)tics-(float)last_tics)/CLOCKS_PER_SEC;
      if (cputime < 0.0) cputime += 4.2949673e+09/(float)CLOCKS_PER_SEC;
      totcputime += cputime;
      avgcputime = totcputime/(step+1);
      //fprintf(stdout,"  took %g sec\n",cputime);
      //fprintf(stdout,"  total %g and mean %g sec\n",totcputime,avgcputime);
      //fprintf(stdout,"mean step time: %.3g sec, effective Re %g\n",
      //   avgcputime,effective_re);

      // write statistics
      if (writeOutput) {
         vmax = compute_and_write_stats(silent,step,dt,simtime,cputime,nx,ny,xbdry,ybdry,u);
      }

      if (first_time) first_time = FALSE;

      // ----------------------------
      step++;
      simtime += dt;
   }

   // close out the stats file properly
   if (writeOutput && maxstep > 0) {
      compute_and_write_stats(silent,-1,dt,simtime,cputime,nx,ny,xbdry,ybdry,u);
   }

   if (!silent) fprintf(stderr,"\nDone.\n");
   exit(0);
}


/*
 * Do just what it says
 */
float compute_and_write_stats(int silent, int step, float dt, float simtime, float cputime,
      int nx, int ny, int xbdry, int ybdry, float ***u) {

   int i,j;
   float ke,base_mult,multx,multy,velsq,vmax,cn;
   static char outfile[80];
   static int initialized = FALSE;
   static FILE *outp;

   // if last time through, close file pointer and RETURN
   if (step == -1) {
      fclose(outp);
      return(0.0);
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
      initialized = TRUE;
   }

   // find kinetic energy
   ke = 0.;
   base_mult = 1.0 / ((nx-1) * (ny-1));
   for (i=0; i<nx; i++) {
      if (i==0 || i==nx-1) multx = base_mult*0.5;
         else multx = base_mult;
      for (j=0; j<ny; j++) {
         if (j==0 || j==ny-1) multy = multx*0.5;
            else multy = multx;
         ke += multy*(u[XV][i][j]*u[XV][i][j] + u[YV][i][j]*u[YV][i][j]);
      }
   }
   ke /= 2.;
   //fprintf(stdout,"ke is %g\n",ke);

   // find vmax
   vmax = 0.;
   for (i=0; i<nx; i++) {
      for (j=0; j<ny; j++) {
         velsq = u[XV][i][j]*u[XV][i][j] + u[YV][i][j]*u[YV][i][j];
         if (velsq > vmax) vmax = velsq;
      }
   }
   vmax = sqrt(vmax);
   cn = vmax*(nx+1)*dt;

   // write the stats here
   if (!silent) {
      fprintf(stdout,"%d %g %g %g %g %g\n",step,simtime,ke,vmax,cn,cputime);
      fflush(stdout);
   }
   fprintf(outp,"%d %g %g %g %g %g\n",step,simtime,ke,vmax,cn,cputime);
   fflush(outp);

   return(vmax);
}


/*
 * Drop a groovy paint splat
 */
void paint_splat (float sx,float sy, float ex,float ey, float rad,float thresh,
                  float r,float g,float b, int nspots,
                  int nx,int ny, float **red,float **grn,float **blu) {

   int i,ix,iy;
   float cx,cy,px,py,temp;
   float **spot;

   // must use isosurface of Gaussian-splattered spots
   // first, create x,y,z locations and strengths of blobs
   spot = allocate_2d_array_f(nspots,3);
   for (i=0; i<nspots; i++) {
      cx = i/(float)nspots;
      // create two Gaussian random numbers for the spot location
      px = (float)(rand())/(float)(RAND_MAX);
      py = (float)(rand())/(float)(RAND_MAX);
      spot[i][0] = sx+(ex-sx)*cx+rad*(sqrt(-2.*log(px))*cos(2.*M_PI*py));
      spot[i][1] = sy+(ey-sy)*cx+rad*(sqrt(-2.*log(px))*sin(2.*M_PI*py));
      // and select a strength for the spot
      py = (float)(rand())/(float)(RAND_MAX);
      spot[i][2] = (cx*(1.-cx))*py/sqrt(-2.*log(px));
   }
   // now, loop over all pixels, perform a summation to see if the
   //   pixel is "inside" or "outside" of the isosurface
   for (ix=0; ix<nx; ix++) {
      for (iy=0; iy<ny; iy++) {
         // pixel position
         px = ix/(float)nx;
         py = iy/(float)ny;
         // scalar potential
         temp = 0.;
         for (i=0; i<nspots; i++) {
            cx = spot[i][0] - px;
            cy = spot[i][1] - py;
            temp += spot[i][2]/(cx*cx+cy*cy);
         }
         //fprintf(stderr,"%d %d  %g\n",ix,iy,temp);
         if (temp > thresh) {
            // this is how to do inky color
            //red[ix][iy] = 1.-4*(1.-r);
            //grn[ix][iy] = 1.-4*(1.-g);
            //blu[ix][iy] = 1.-4*(1.-b);
            // do the color sent
            red[ix][iy] = r;
            grn[ix][iy] = g;
            blu[ix][iy] = b;
         }
      }
   }
}


/*
 * grab a color from the color image
 */
int get_random_color (float ***c, int nx, int ny, float *thisc) {

   int ix,iy;

   // choose random coordinate
   ix = nx*(float)rand()/(float)RAND_MAX;
   iy = ny*(float)rand()/(float)RAND_MAX;

   thisc[0] = c[0][ix][iy];
   thisc[1] = c[1][ix][iy];
   thisc[2] = c[2][ix][iy];
}


/*
 * This function writes basic usage information to stderr,
 * and then quits. Too bad.
 */
int Usage(char progname[80],int status) {

   static char **cpp, *help_message[] = {
   "where [-options] are one or more of the following:                         ",
   "                                                                           ",
   "   -x [int]    resolution in x-direction; default=512                      ",
   "   -y [int]    resolution in y-direction; default=512                      ",
   "                                                                           ",
   "   -px         make domain periodic in x-direction; default is wall bc     ",
   "   -py         make domain periodic in y-direction; default is wall bc     ",
   "   -open       make domain open in both directions; default is wall bc     ",
   "   -bcl [float] set wall-tangent bondary vel on left side                  ",
   "   -bcr [float] set wall-tangent bondary vel on right side                 ",
   "   -bcb [float] set wall-tangent bondary vel on bottom side                ",
   "   -bct [float] set wall-tangent bondary vel on top side                   ",
   "                                                                           ",
   "   -fs [float] [float] set freestream in x and y directions, sets open bc  ",
   "                                                                           ",
   "   -uprint     write velocity field (vel_00000.png RGB )                   ",
   "   -uscale [float] scale factor for velocity output; default is 1.0        ",
   "                                                                           ",
   "   -vprint     write vorticity field (vort_00000.png)                      ",
   "   -vscale [float] scale factor for vorticity output; default is 10.       ",
   "                                                                           ",
   "   -vf [name]  read PNG file and use as initial vorticity field            ",
   "   -vd [float] set momentum diffusivity; default is 1.e-3                  ",
   "                                                                           ",
   "   -vdf [name] [float] [float] [float]                                     ",
   "               read PNG file and set variable momentum diffusivity, floats ",
   "               are min, max diffusivity, viscosity diffusivity             ",
   "   -muprint    write momentum viscosity field (mu_00000.png)               ",
   "                                                                           ",
   "   -t          track temperature/density field                             ",
   "   -tprint     write temperature/density field (temp_00000.png)            ",
   "   -tf [name]  read PNG file and use as temperature/density field          ",
   "   -td [float] set temperature/density diffusivity; default is 1.e-3       ",
   "   -b [float]  set the Boussinesq number; default=1, turns off -dr         ",
   "   -dr [float] set the absolute density ratio; def=10., turns off -b       ",
   "                                                                           ",
   "   -qf [name]  read PNG file and use as heat source field                  ",
   "   -qc [float] set the scaling on the heat source term                     ",
   "                                                                           ",
   "   -c          track and print color scalars (out_00000.png)               ",
   "   -cf [name]  read PNG file and use as color scalar field                 ",
   "   -cd [float] set color diffusivity, absolute scale; default is 1.e-3     ",
   "   -cdf [name] read PNG file and set color diffusivity                     ",
   "                                                                           ",
   "   -div        track and print dilation field (div_00000.png)              ",
   "   -df [name]  read PNG file and use as initial dilation/divergence field  ",
   "                                                                           ",
   "   -m          track mask field                                            ",
   "   -mf [name]  read PNG file and use as flow mask, black areas are solid,  ",
   "               white areas are completely open                             ",
   "   -me [float] set residual error in iterations; default is 1.e-3          ",
   "                                                                           ",
   "   -dt [float] set the time step increment, if unused, Courant number will ",
   "               determine the time step; no default                         ",
   "   -cn [float] set the Courant number, based on a speed of 1.0; try        ",
   "               10 to 50 for a small simulation (32 to 128) or 50 to 500    ",
   "               for a larger one (256 to 2048); default=10                  ",
   "                                                                           ",
   "   -ores [int] set output resolution to num times simulation resolution;   ",
   "               default=1                                                   ",
   "                                                                           ",
   "   -every [int] write only every so many steps, default=1 (every step)     ",
   "   -noout      suppress all file output                                    ",
   "   -step [int] run a fixed number of steps; default=99999                  ",
   "                                                                           ",
   "   -stam       use Stam's velocity method-of-characteristics instead of    ",
   "               the default vorticity m.o.c.                                ",
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

