/*
 * vic2d.c - driver file for two-dimensional vortex-in-cell method-of-characteristics
 *
 * Copyright 2004-23 Mark J. Stock <mstock@umich.edu>
 *
 * a 2D vortex-in-cell method which uses the method of characteristics for the
 * convection step, and a single explicit step for diffusion and vorticity
 * creation from walls and masks
 *
 */

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
#include "maskops.h"
#include "inout.h"
#include "particles.h"
#include "vicmoc.h"

float compute_and_write_stats(int, int, float, float, double, int, int, int, int, float***,float**);
void paint_splat (float,float,float,float,float,float,float,float,float,int,int,int,float**,float**,float**);
void get_random_color (float***,int,int,float*);
void get_color (float***,int,int,float,float,float*);
int Usage(char[MAXCHARS],int);

int main(int argc,char **argv) {

   static bool first_time = true;
   bool isStam = false;
   int nx = 513;
   int ny = 513;
   int cnx,cny,cmn = 0;
   int xbdry = WALL;
   int ybdry = WALL;
   int step = -1;
   int maxstep = 99999;
   int writeevery = 1;
   int checkpointevery = -1;
   int outscale = 1;
   bool use_16bpp = false;
   bool silent = false;

   float sc_diffus[MAX_SCALARS];
   bool use_TEMP = false;
   bool use_DF = false;
   bool use_COLOR = false;
   bool reread_color = false;
   bool use_cint = false;
   int cint_interval = 100;
   bool darkenonly = true;
   bool freeze_flow = false;
   int freeze_at = -1;
   bool stop_flow = false;
   int stop_at = -1;
   bool use_MASK = false;
   bool use_image = false;
   bool use_strong_strat = false;
   bool overlay_color_lr = false;
   bool overlay_color_tb = false;
   bool vort_reread = false;
   bool vort_suppress = false;
   bool adaptive_mask = false;

   bool writeOutput = true;
   bool print_vort = false;
   bool print_vel = false;
   bool print_temp = false;
   bool print_mu = false;
   bool print_mask = false;

   float md = 1.e-3;			// dimensional momentum diffusivity
   float mdlow = -1.0;
   float mdhigh = -1.0;
   float vd = 1.e-3;			// dimensional viscosity diffusivity
   float td = 1.e-3;			// dimensional density diffusivity
   float tdlow = -1.0;
   float tdhigh = -1.0;
   float cd = 1.e-3;			// dimensional color diffusivity
   float cdlow = -1.0;
   float cdhigh = -1.0;
   float bn = 1.0;				// Boussinesq coefficient
   float dens_ratio = 1.0;		// ratio of max to min density
   float heat_coeff = 1.0;		// coefficient on the heat term
   float dt = -1.0;				// time step size
   float courant = 10.;			// fixed non-dim time step size
   float courantconst = -1.;	// target non-dim time step size
   float freestream[2] = {0.0, 0.0};
   bool dynamic_freestream = false;
   float fs_start[2] = {0.0, 0.0};
   int fs_start_step = 0;
   float fs_end[2] = {0.0, 0.0};
   int fs_end_step = 100000;
   float wallvel[4] = {0., 0., 0., 0.};	// left, right, bottom, top
   bool recalc_vel = true;
   bool move_colors = true;
   int ccnidx = 0;
   int gravtype = 0;
   float gravity[2] = {0.0, 1.0};
   float yf;
   float effective_re;
   float px,py,cx,cy,rad,temp,temp2,maxval,minval,scale;
   float vortscale = 10.;
   float velscale = 1.0;
   float ccnvmax[VMAXAVG];
   float randvortscale = -1;
   float vort_suppress_factor = -1;
   float overlay_fraction = 0.01;
   double walltime = 0.;
   float simtime = 0.;
   float thisc[3];
   float maskerr = 5.e-3;
   float **color_left = NULL;
   float **color_right = NULL;
   float **color_top = NULL;
   float **color_bottom = NULL;

   bool vort_add_rand = false;
   float vort_add_factor = -1;
   int vort_add_start = -1;
   int vort_add_peak = -1;
   int vort_add_down = -1;
   int vort_add_end = -1;

   bool use_PARTICLES = false;
   int pnx = -1;
   int pny = -1;
   int part_draw_fade_frames = 100;
   int part_add_step_start = 0;
   int part_add_step_end = 0;
   int part_add_count_end = 10000;
   int part_rem_step_start = 999999;
   int part_rem_step_end = 999999;
   int part_rem_count_end = 0;
   float particle_speed = 1.0;
   float part_ballistic = 1.0;
   float part_draw_fac = 1.0;
   float part_draw_mass_pow = 0.0;
   float part_draw_vel_pow = 1.0;
   float part_img_diffus = -1.0;

   bool dynamic_mask = false;
   int dmask_step_start = 0;
   int dmask_step_end = 1000;
   float dmask_val_start = 0.5;
   float dmask_val_end = 0.5;
   float dmask_width = 0.5;
   float dmask_power = 1.0;
   float mask_color_mult = 1.0;

   float **u[2];			// velocities
   float **a[MAX_SCALARS];		// other scalars
   float **t[MAX_SCALARS];		// temporary velocities, vorticity, etc
   float **mask;			// flow mask
   float **heat;			// constant heat source map
   float **c[3];			// color image storage for temp field
   float **cint[3];			// integrated color image storage
   float **pc[3];			// particle color image
   float *cml[3];			// linear color map storage
   float **cmi[3];			// 2D color input storage
   float **acc[2];			// Lagrangian acceleration
   float **shear;			// shear magnitude

   struct Particles pts;		// collection of particles

   // bookkeeping
   unsigned long int tics,last_tics;
   struct timeval t_curr, t_last;
   char progname[MAXCHARS];
   char outfileroot[MAXCHARS];
   bool use_vort_img = false;
   char vortfilename[MAXCHARS];
   bool use_temp_img = false;
   char tempfilename[MAXCHARS];
   bool use_heat_img = false;
   char heatfilename[MAXCHARS];
   bool use_color_img = false;
   char colorfilename[MAXCHARS];
   bool use_div_img = false;
   char divfilename[MAXCHARS];
   bool use_mask_img = false;
   char maskfilename[MAXCHARS];
   bool use_color_linear = false;
   bool use_color_area = false;
   char colorsrcfilename[MAXCHARS];
   char cintfileroot[MAXCHARS];
   char mdfilename[MAXCHARS];
   char tdfilename[MAXCHARS];
   char cdfilename[MAXCHARS];

   // Initialize arrays and particles
   for (int i=0; i<MAX_SCALARS; i++) sc_diffus[i] = 1.e-3;
   for (int i=0; i<MAX_SCALARS; i++) a[i] = NULL;
   for (int i=0; i<MAX_SCALARS; i++) t[i] = NULL;
   (void) init_particles(&pts, 10000);

   // read command-line
   (void) strcpy(progname,argv[0]);
   if (argc < 1) (void) Usage(progname,0);
   for (int i=1; i<argc; i++) {

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
         use_vort_img = true;
      } else if (strncmp(argv[i], "-vrr", 4) == 0) {
         vort_reread = true;
      } else if (strncmp(argv[i], "-vdf", 4) == 0) {
         strcpy (mdfilename,argv[++i]);
         mdlow = atof(argv[++i]);
         mdhigh = atof(argv[++i]);
         vd = atof(argv[++i]);
      } else if (strncmp(argv[i], "-muprint", 3) == 0) {
         print_mu = true;
      } else if (strncmp(argv[i], "-vd", 3) == 0) {
         md = atof(argv[++i]);
      } else if (strncmp(argv[i], "-vas", 4) == 0) {
         vort_add_start = atoi(argv[++i]);
         vort_add_peak = atoi(argv[++i]);
         vort_add_down = atoi(argv[++i]);
         vort_add_end = atoi(argv[++i]);
      } else if (strncmp(argv[i], "-va", 3) == 0) {
         vort_add_rand = true;
         vort_add_factor = atof(argv[++i]);
      } else if (strncmp(argv[i], "-vr", 3) == 0) {
         vort_suppress = true;
         vort_suppress_factor = atof(argv[++i]);
      } else if (strncmp(argv[i], "-vprint", 3) == 0) {
         print_vort = true;
      } else if (strncmp(argv[i], "-vscale", 3) == 0) {
         vortscale = atof(argv[++i]);

      } else if (strncmp(argv[i], "-uscale", 3) == 0) {
         velscale = atof(argv[++i]);
      } else if (strncmp(argv[i], "-uprint", 3) == 0) {
         print_vel = true;
      } else if (strncmp(argv[i], "-freeze", 3) == 0) {
         freeze_flow = true;
         if (argc > i+1) {
            if (isdigit((int)argv[i+1][0]) || isdigit((int)argv[i+1][1])) {
               // first number after keyword is frame at which to freeze
               freeze_at = atoi(argv[++i]);
            }
         }
      } else if (strncmp(argv[i], "-stop", 4) == 0) {
         stop_flow = true;
         if (argc > i+1) {
            if (isdigit((int)argv[i+1][0]) || isdigit((int)argv[i+1][1])) {
               // first number after keyword is frame at which to stop
               stop_at = atoi(argv[++i]);
            }
         }
      } else if (strncmp(argv[i], "-constcn", 8) == 0) {
         courantconst = atof(argv[++i]);
      } else if (strncmp(argv[i], "-checkpoint", 4) == 0) {
         checkpointevery = atoi(argv[++i]);
         if (checkpointevery < 1) checkpointevery = -1;

      } else if (strncmp(argv[i], "-tf", 3) == 0) {
         strcpy (tempfilename,argv[++i]);
         use_temp_img = true;
         use_TEMP = true;
      } else if (strncmp(argv[i], "-qf", 3) == 0) {
         strcpy (heatfilename,argv[++i]);
         use_heat_img = true;
         use_TEMP = true;
      } else if (strncmp(argv[i], "-qc", 3) == 0) {
         heat_coeff = atof(argv[++i]);
         use_TEMP = true;
      } else if (strncmp(argv[i], "-tdf", 4) == 0) {
         strcpy (tdfilename,argv[++i]);
         tdlow = atof(argv[++i]);
         tdhigh = atof(argv[++i]);
      } else if (strncmp(argv[i], "-td", 3) == 0) {
         td = atof(argv[++i]);
      } else if (strncmp(argv[i], "-tprint", 3) == 0) {
         print_temp = true;
         use_TEMP = true;
      } else if (strncmp(argv[i], "-t", 2) == 0) {
         use_TEMP = true;
      } else if (strncmp(argv[i], "-b", 2) == 0) {
         bn = atof(argv[++i]);
         use_strong_strat = false;
      } else if (strncmp(argv[i], "-dr", 3) == 0) {
         dens_ratio = atof(argv[++i]);
         use_strong_strat = true;

      } else if (strncmp(argv[i], "-cf", 3) == 0) {
         strcpy (colorfilename,argv[++i]);
         use_color_img = true;
         use_COLOR = true;
      } else if (strncmp(argv[i], "-cdf", 4) == 0) {
         strcpy (cdfilename,argv[++i]);
         cdlow = atof(argv[++i]);
         cdhigh = atof(argv[++i]);
      } else if (strncmp(argv[i], "-cr", 3) == 0) {
         reread_color = true;
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
               //   darkenonly = true;
               //   overlay_fraction *= -1.;
               //}
            }
         }
      } else if (strncmp(argv[i], "-cd", 3) == 0) {
         cd = atof(argv[++i]);
      } else if (strncmp(argv[i], "-cl", 3) == 0) {
         strcpy (colorsrcfilename,argv[++i]);
         use_color_linear = true;
      } else if (strncmp(argv[i], "-ca", 3) == 0) {
         strcpy (colorsrcfilename,argv[++i]);
         use_color_area = true;
      } else if (strncmp(argv[i], "-cint", 5) == 0) {
         use_cint = true;
         cint_interval = atoi(argv[++i]);
         if (cint_interval < 2) cint_interval = 2;
      } else if (strncmp(argv[i], "-c", 2) == 0) {
         use_COLOR = true;

      } else if (strncmp(argv[i], "-df", 3) == 0) {
         strcpy (divfilename,argv[++i]);
         use_div_img = true;
      } else if (strncmp(argv[i], "-div", 3) == 0) {
         use_DF = true;

      } else if (strncmp(argv[i], "-mf", 3) == 0) {
         strcpy (maskfilename,argv[++i]);
         use_mask_img = true;
         use_MASK = true;
      } else if (strncmp(argv[i], "-me", 3) == 0) {
         maskerr = atof(argv[++i]);
      } else if (strncmp(argv[i], "-mprint", 4) == 0) {
         print_mask = true;
      } else if (strncmp(argv[i], "-mdyn", 3) == 0) {
         use_MASK = true;
         dynamic_mask = true;
         strcpy (maskfilename,argv[++i]);
         dmask_step_start = atoi(argv[++i]);
         dmask_val_start = atof(argv[++i]);
         dmask_step_end = atoi(argv[++i]);
         dmask_val_end = atof(argv[++i]);
         dmask_width = atof(argv[++i]);
      } else if (strncmp(argv[i], "-mpow", 4) == 0) {
         dmask_power = atof(argv[++i]);
      } else if (strncmp(argv[i], "-madaptu", 3) == 0) {
         adaptive_mask = true;
      } else if (strncmp(argv[i], "-mcm", 3) == 0) {
         mask_color_mult = atof(argv[++i]);
      } else if (strncmp(argv[i], "-m", 2) == 0) {
         use_MASK = true;

      } else if (strncmp(argv[i], "-pa", 3) == 0) {
         use_PARTICLES = true;
         part_add_step_start = atoi(argv[++i]);
         part_add_step_end = atoi(argv[++i]);
         part_add_count_end = atoi(argv[++i]);
      } else if (strncmp(argv[i], "-pr", 3) == 0) {
         part_rem_step_start = atoi(argv[++i]);
         part_rem_step_end = atoi(argv[++i]);
         part_rem_count_end = atoi(argv[++i]);
      } else if (strncmp(argv[i], "-ps", 3) == 0) {
         particle_speed = atof(argv[++i]);
      } else if (strncmp(argv[i], "-pb", 3) == 0) {
         part_ballistic = atof(argv[++i]);
      } else if (strncmp(argv[i], "-pf", 3) == 0) {
         part_draw_fade_frames = atof(argv[++i]);
      } else if (strncmp(argv[i], "-pd", 3) == 0) {
         part_draw_fac = atof(argv[++i]);
         part_draw_mass_pow = atof(argv[++i]);
         part_draw_vel_pow = atof(argv[++i]);
      } else if (strncmp(argv[i], "-pid", 4) == 0) {
         part_img_diffus = atof(argv[++i]);
      } else if (strncmp(argv[i], "-pnx", 4) == 0) {
         pnx = atoi(argv[++i]);
      } else if (strncmp(argv[i], "-pny", 4) == 0) {
         pny = atoi(argv[++i]);

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
         isStam = true;
      } else if (strncmp(argv[i], "-fsdyn", 4) == 0) {
         dynamic_freestream = true;
         fs_start_step = atoi(argv[++i]);
         fs_start[0] = atof(argv[++i]);
         fs_start[1] = atof(argv[++i]);
         fs_end_step = atoi(argv[++i]);
         fs_end[0] = atof(argv[++i]);
         fs_end[1] = atof(argv[++i]);
      } else if (strncmp(argv[i], "-fs", 3) == 0) {
         freestream[0] = atof(argv[++i]);
         freestream[1] = atof(argv[++i]);
      } else if (strncmp(argv[i], "-grav", 2) == 0) {
         gravtype = 0;
         gravity[0] = atof(argv[++i]);
         gravity[1] = atof(argv[++i]);
      } else if (strncmp(argv[i], "-rgrav", 6) == 0) {
         gravtype = 1;
         gravity[0] = atof(argv[++i]);
         gravity[1] = atof(argv[++i]);


      } else if (strncmp(argv[i], "-randvortscale", 6) == 0) {
         randvortscale = atof(argv[++i]);

      } else if (strncmp(argv[i], "-noout", 6) == 0) {
         writeOutput = false;
      } else if (strncmp(argv[i], "-q", 2) == 0) {
         silent = true;
      } else if (strncmp(argv[i], "-8", 2) == 0) {
         use_16bpp = false;
      } else if (strncmp(argv[i], "-16", 3) == 0) {
         use_16bpp = true;
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

   // save integrated color field?
   if (use_cint && use_COLOR) {
      cint[0] = allocate_2d_array_f(nx,ny);
      cint[1] = allocate_2d_array_f(nx,ny);
      cint[2] = allocate_2d_array_f(nx,ny);
      for (int ix=0; ix<nx; ix++) { for (int iy=0; iy<ny; iy++) { cint[0][ix][iy] = 0.0f; } }
      for (int ix=0; ix<nx; ix++) { for (int iy=0; iy<ny; iy++) { cint[1][ix][iy] = 0.0f; } }
      for (int ix=0; ix<nx; ix++) { for (int iy=0; iy<ny; iy++) { cint[2][ix][iy] = 0.0f; } }
   } else {
      use_cint = false;
   }

   // variable viscosities

   // variable momentum viscosity
   if (mdlow > 0.0 && mdhigh > 0.0) {
      if (!silent) fprintf(stdout,"Using variable momentum viscosity\n");
      // viscosity diffuses with its own special (constant) diffusivity,
      sc_diffus[MD] = vd;
      a[MD] = allocate_2d_array_f(nx,ny);
      // read grayscale PNG of exactly nx by ny resolution
      read_png(mdfilename,nx,ny,false,false,1.0,false,
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
   for (int i=0; i<MAX_SCALARS; i++) {
      if (a[i] != NULL) {
         t[i] = allocate_2d_array_f(nx,ny);
         if (!silent) fprintf(stdout,"Array %d has diffusivity %g\n",i,sc_diffus[i]);
      }
   }

   // -----------------------------------------------
   // Set initial conditions

   // Set vorticity ----------------------------------

   // first, zero it
   for (int ix=0; ix<nx; ix++) {
      for (int iy=0; iy<ny; iy++) {
         a[W2][ix][iy] = 0.0;
      }
   }

   if (false) {
      // create a wavy shear layer of positive vorticity
      for (int ix=0; ix<nx; ix++) {
         px = (float)ix/(float)(nx-1);
         for (int iy=0; iy<ny; iy++) {
            py = (float)iy/(float)(nx-1);
            a[W2][ix][iy] = 10.0 - 20.0*fabs((yf/2.0) - py - 0.05*sin(px*2.0*M_PI));
            if (a[W2][ix][iy] < 0.0) a[W2][ix][iy] = 0.0;
         }
      }
   } else if (false) {
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
   } else if (false) {
      // create a random field of vorticity
      //scale = 20.;
      scale = 100.;
      for (int ix=0; ix<nx; ix++) {
         for (int iy=0; iy<ny; iy++) {
            a[W2][ix][iy] = scale*rand()/RAND_MAX - scale*0.5;
         }
      }
   } else if (false) {
      add_smooth_circular_blob(nx,ny,xbdry,ybdry,a[W2],0.6,0.2,0.1,10.0);
      add_smooth_circular_blob(nx,ny,xbdry,ybdry,a[W2],0.4,0.16,0.1,-10.0);
   } else if (false) {
      add_smooth_circular_blob(nx,ny,xbdry,ybdry,a[W2],0.57,0.37,0.05,50.0);
      add_smooth_circular_blob(nx,ny,xbdry,ybdry,a[W2],0.41,0.34,0.05,-50.0);
   } else if (false) {
      // make the vorticity defined in Minion and Brown, JCP 138
      for (int ix=0; ix<nx; ix++) {
         px = (float)ix/(float)(nx-1);
         for (int iy=0; iy<ny; iy++) {
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
   } else if (false) {
      // make the vorticity defined in Koumoutsakos, JCP, 1997, type I
      for (int ix=0; ix<nx; ix++) {
         px = (float)ix/(float)(nx-1);
         for (int iy=0; iy<ny; iy++) {
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
   } else if (false) {
      // create a sinusoidal interface
      //float thickness = 0.001;
      float thickness = 2.0/(float)nx;
      for (int ix=0; ix<nx; ix++) {
         px = (ix - (float)nx/2.)/(float)nx;
         yf = 0.02*sin(px*9.*M_PI);
         for (int iy=0; iy<ny; iy++) {
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
      // inside read_png we set 127 or 32767 to be 0.0
      read_png(vortfilename,nx,ny,false,
               false,1.0,false,
               a[W2],-vortscale,2.0*vortscale,
               NULL, -vortscale,2.0*vortscale,
               NULL, -vortscale,2.0*vortscale);
   }
   if (randvortscale > 0.0) {
      // create a random field of vorticity
      for (int ix=0; ix<nx; ix++) {
         for (int iy=0; iy<ny; iy++) {
            a[W2][ix][iy] += 2.0*randvortscale*(rand()/(float)RAND_MAX - 0.5);
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
      for (int ix=0; ix<nx; ix++) {
         for (int iy=0; iy<ny; iy++) {
            // use 0.0 for all solid
            mask[ix][iy] = 1.0;
         }
      }
      if (false) {
         // create a circular mask at cx,cy
         cx = 0.3;
         cy = 0.5;
         rad = 0.15;
         //temp2 = 0.0025;
         temp2 = 2.0/(float)nx;
         rad -= temp2/2.;
         for (int ix=0; ix<nx; ix++) {
            px = (float)ix/(float)(nx-1);
            for (int iy=0; iy<ny; iy++) {
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
      if (false) {
         // create a solar tower
         cx = 0.5;	// x center
         cy = 0.0;	// y center
         minval = 0.01; // inner surface
         maxval = 0.06; // outer surface
         temp2 = 0.002;	// thickness
         scale = 0.3;	// length
         for (int ix=0; ix<nx; ix++) {
           px = fabs((float)ix/(float)(nx-1) - cx);
           if (px < scale) {
            for (int iy=0; iy<ny; iy++) {
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
         read_png(maskfilename,nx,ny,false,false,1.0,false,
            mask,0.0,1.0,NULL,0.0,1.0,NULL,0.0,1.0);
      }
      // normalize mask, and set 1.0=solid, 0.0=open (input is opposite)
      maxval = -9.9e+9;
      minval = 9.9e+9;
      for (int ix=0; ix<nx; ix++) {
         for (int iy=0; iy<ny; iy++) {
            if (mask[ix][iy] > maxval) maxval = mask[ix][iy];
            if (mask[ix][iy] < minval) minval = mask[ix][iy];
         }
      }
      temp = maxval-minval;
      if (temp > 1.e-5) {
      for (int ix=0; ix<nx; ix++) {
         for (int iy=0; iy<ny; iy++) {
            mask[ix][iy] = (maxval-mask[ix][iy])/temp;
            if (mask[ix][iy] < 0.0) mask[ix][iy] = 0.0;
            if (mask[ix][iy] > 1.0) mask[ix][iy] = 1.0;
         }
      }
      }

      if (dynamic_mask) {
         (void) set_mask_from_temporal (0, nx, ny, mask, maskfilename,
                                        dmask_step_start, dmask_val_start,
                                        dmask_step_end, dmask_val_end,
                                        dmask_width, dmask_power);
      }
   }

   // optionally generate repeatedly-overlaid mask
   //(void) overlay_mask (nx, ny, mask);

   // set freestream ----------------------------------

   if (dynamic_freestream) {
      if (0 < fs_start_step) {
         freestream[0] = fs_start[0];
         freestream[1] = fs_start[1];
      } else if (0 < fs_end_step) {
         const float factor = (float)(0-fs_start_step) / (float)(fs_end_step-fs_start_step);
         freestream[0] = fs_start[0] + factor*(fs_end[0]-fs_start[0]);
         freestream[1] = fs_start[1] + factor*(fs_end[1]-fs_start[1]);
      } else {
         freestream[0] = fs_end[0];
         freestream[1] = fs_end[1];
      }
   }

   // Initial velocity solve -------------------------

   // if you want to show velocity on the first step, solve for it here
   find_vels_2d (silent,step,isStam,nx,ny,xbdry,ybdry,freestream,wallvel,u[XV],u[YV],a[W2],use_MASK,mask,maskerr);

   // find vmax (might need this)
   float vmax = 0.;
   for (int i=0; i<nx; i++) {
      for (int j=0; j<ny; j++) {
         float velsq = u[XV][i][j]*u[XV][i][j] + u[YV][i][j]*u[YV][i][j];
         if (velsq > vmax) vmax = velsq;
      }
   }
   vmax = sqrt(vmax);
   // initialize the vmax averaging array
   vmax = courantconst / (dt*(nx+1));
   for (int i=0; i<VMAXAVG; i++) ccnvmax[i] = vmax;


   // set color ---------------------------------------

   // initialize an array of blocks to add during the run
   if (false) {
      populate_block_array(nx,ny);
   }

   // first, find what index the scalar uses
   if (use_COLOR) {
      // zero the arrays
      for (int ix=0; ix<nx; ix++) for (int iy=0; iy<ny; iy++) a[RR][ix][iy] = 0.0;
      for (int ix=0; ix<nx; ix++) for (int iy=0; iy<ny; iy++) a[GG][ix][iy] = 0.0;
      for (int ix=0; ix<nx; ix++) for (int iy=0; iy<ny; iy++) a[BB][ix][iy] = 0.0;
      //for (int ix=0; ix<nx; ix++) for (int iy=0; iy<ny; iy++) a[RR][ix][iy] = 1.0;
      //for (int ix=0; ix<nx; ix++) for (int iy=0; iy<ny; iy++) a[GG][ix][iy] = 1.0;
      //for (int ix=0; ix<nx; ix++) for (int iy=0; iy<ny; iy++) a[BB][ix][iy] = 1.0;
      //for (ix=0; ix<nx; ix++) for (int iy=0; iy<ny; iy++) a[RR][ix][iy] = 0.10546875;
      //for (ix=0; ix<nx; ix++) for (int iy=0; iy<ny; iy++) a[GG][ix][iy] = 0.15234375;
      //for (ix=0; ix<nx; ix++) for (int iy=0; iy<ny; iy++) a[BB][ix][iy] = 0.05859375;
      if (use_color_img) {
         // read color PNG of exactly nx by ny resolution into
         //   the given fields
         read_png (colorfilename,nx,ny,true,false,1.0,false,
            a[RR],0.0,1.0,a[GG],0.0,1.0,a[BB],0.0,1.0);
         //read_png (colorfilename,nx,ny,true,false,1.0,false,
         //   a[RR],-3.0,4.0,a[GG],-3.0,4.0,a[BB],-3.0,4.0);
      }
      if (false) {
         // left is red, right is green
         for (int ix=0; ix<nx/2; ix++) {
            for (int iy=0; iy<ny; iy++) {
               a[RR][ix][iy] = 1.0;
               a[GG][ix][iy] = 0.91372549;
               a[BB][ix][iy] = 0.168627451;
            }
         }
         for (int ix=nx/2; ix<nx; ix++) {
            for (int iy=0; iy<ny; iy++) {
               a[RR][ix][iy] = 0.1171875;
               a[GG][ix][iy] = 0.0546875;
               a[BB][ix][iy] = 0.7265625;
            }
         }
      }
      if (false) {
         // first, force background to one color
         for (int ix=0; ix<nx; ix++) for (int iy=0; iy<ny; iy++) a[RR][ix][iy] = 235./256.;
         for (int ix=0; ix<nx; ix++) for (int iy=0; iy<ny; iy++) a[GG][ix][iy] = 245./256.;
         for (int ix=0; ix<nx; ix++) for (int iy=0; iy<ny; iy++) a[BB][ix][iy] = 253./256.;
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
      if (false) {
         // scale input vorticity by pixel brightness
         for (int ix=0; ix<nx; ix++) {
            for (int iy=0; iy<ny; iy++) {
               float brite = 0.3*a[RR][ix][iy] + 0.6*a[GG][ix][iy] + 0.1*a[BB][ix][iy];
               //a[W2][ix][iy] *= 1.0-brite;
               a[W2][ix][iy] *= MAX(0.0, brite-0.1);
            }
         }
      }

      // save the edge colors
      for (int ix=0; ix<nx; ix++) {
         color_bottom[0][ix] = a[RR][ix][0];
         color_bottom[1][ix] = a[GG][ix][0];
         color_bottom[2][ix] = a[BB][ix][0];
         color_top[0][ix] = a[RR][ix][ny-1];
         color_top[1][ix] = a[GG][ix][ny-1];
         color_top[2][ix] = a[BB][ix][ny-1];
      }
      for (int iy=0; iy<ny; iy++) {
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
      for (int ix=0; ix<nx; ix++) {
         for (int iy=0; iy<ny; iy++) {
            a[SF][ix][iy] = 0.0;
         }
      }

      if (use_temp_img) {
         // read grayscale PNG of exactly nx by ny resolution
         //read_png(tempfilename,nx,ny,false,true,1.0,false,
         read_png(tempfilename,nx,ny,false,false,1.0,false,
            a[SF],-1.0,2.0,NULL,0.0,1.0,NULL,0.0,1.0);
      }

      if (use_heat_img) {
         // read grayscale PNG of exactly nx by ny resolution
         read_png(heatfilename,nx,ny,false,false,1.0,false,
            heat,-1.0,2.0,NULL,0.0,1.0,NULL,0.0,1.0);

         for (int ix=0; ix<nx; ix++) {
            for (int iy=0; iy<ny; iy++) {
               a[SF][ix][iy] += heat[ix][iy];
               if (a[SF][ix][iy] > 1.0) a[SF][ix][iy] = 1.0;
               if (a[SF][ix][iy] < -1.0) a[SF][ix][iy] = -1.0;
            }
         }
      }

      if (false) {
         // create a sharp blob of scalar
         for (int ix=0; ix<nx; ix++) {
            px = (float)ix/(float)(nx-1);
            for (int iy=0; iy<ny; iy++) {
               py = (float)iy/(float)(nx-1);
               if (sqrt(pow(px-0.4,2)+pow(py-0.15,2)) < 0.1)
                  a[SF][ix][iy] = 10.0;
            }
         }
      }
      if (false) {
         // create a smooth blob of scalar
         // add_smooth_circular_blob(nx,ny,xbdry,ybdry,a[SF],0.7,0.15,0.1,1.0);
         add_smooth_circular_blob(nx,ny,xbdry,ybdry,a[SF],0.4,0.0,0.07,-1.0);
      }
      if (false) {
         // create bands of scalar
         for (int ix=0; ix<nx; ix++) {
            for (int iy=0; iy<ny; iy++) {
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
      if (false) {
         // create square bands of scalar
         for (int ix=0; ix<nx; ix++) {
            for (int iy=0; iy<ny; iy++) {
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
      if (false) {
         // create a RTI with a sinusoidal interface
         //float thickness = 0.001;
         float thickness = 2.0/(float)nx;
         //float thickness = 1.0/(float)nx;
         for (int ix=0; ix<nx; ix++) {
            px = (ix - (float)nx/2.)/(float)nx;
            //yf = 0.02*sin(px*9.*M_PI);
            yf = 0.02*sin(px*8.*M_PI);
            for (int iy=0; iy<ny; iy++) {
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
      if (false) {
         // create another sinusoidal interface (bottom one)
         //float thickness = 0.001;
         float thickness = 2.0/(float)nx;
         for (int ix=0; ix<nx; ix++) {
            px = (ix - (float)nx/2.)/(float)nx;
            yf = 0.02*sin(px*9.*M_PI) - 0.02;
            for (int iy=0; iy<ny; iy++) {
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
      if (false) {
        // create a random field of scalar
        scale = 10.;
        for (int ix=0; ix<nx; ix++) {
          for (int iy=0; iy<ny; iy++) {
            a[SF][ix][iy] += scale*rand()/RAND_MAX - 0.5*scale + 0.5;
          }
        }
      }
      if (false && use_COLOR) {
         // use the brightness component to indicate density!
         for (int ix=0; ix<nx; ix++) {
            for (int iy=0; iy<ny; iy++) {
               a[SF][ix][iy] = 1.0 - 0.6*a[GG][ix][iy] -
                                0.3*a[RR][ix][iy] -
                                0.1*a[BB][ix][iy];
            }
         }
      }
   }

   // Set color map -------------------------------------------

   // load in an image containing colors to grab
   if (use_color_linear || use_color_area) {

      // interrogate the header for resolution
      read_png_res (colorsrcfilename,&cny,&cnx);

      // allocate space and read image
      cmi[0] = allocate_2d_array_f(cnx,cny);
      cmi[1] = allocate_2d_array_f(cnx,cny);
      cmi[2] = allocate_2d_array_f(cnx,cny);
      read_png (colorsrcfilename,cnx,cny,true,false,1.0,false,
                cmi[0],0.0,1.0,cmi[1],0.0,1.0,cmi[2],0.0,1.0);

      // later on, grab colors with
      //(void) get_random_color (cmi,cnx,cny,thisc);
   }

   if (use_color_linear) {
      // set the linear color map (not a 2D image) along the longest side
      cmn = 0;
      if (cnx > cny) {
         cml[0] = allocate_1d_array_f(cnx);
         cml[1] = allocate_1d_array_f(cnx);
         cml[2] = allocate_1d_array_f(cnx);
         float factor = 1.0;
         for (int ix=0; ix<cnx; ix++) {
            //factor = 1.2 - 2.0*pow(0.6 - (float)ix/(float)(cnx-1), 2);
            //factor = 1.4 - 3.0*pow(0.6 - (float)ix/(float)(cnx-1), 2);
            factor = 1.6 - 3.5*pow(0.65 - (float)ix/(float)(cnx-1), 2);
            cml[0][ix] = factor*cmi[0][ix][0];
            cml[1][ix] = factor*cmi[1][ix][0];
            cml[2][ix] = factor*cmi[2][ix][0];
         }
         cmn = cnx;
      } else {
         cml[0] = allocate_1d_array_f(cny);
         cml[1] = allocate_1d_array_f(cny);
         cml[2] = allocate_1d_array_f(cny);
         float factor = 1.0;
         for (int iy=0; iy<cny; iy++) {
            //factor = 1.4 - 3.0*pow(0.6 - (float)iy/(float)(cny-1), 2);
            cml[0][iy] = factor*cmi[0][0][iy];
            cml[1][iy] = factor*cmi[1][0][iy];
            cml[2][iy] = factor*cmi[2][0][iy];
         }
         cmn = cny;
      }
   }

   // Set particles -------------------------------------------

   //if (use_PARTICLES && (use_color_linear || use_color_area)) {
   if (use_PARTICLES) {
      //(void) add_block_of_particles (&pts, 200000, 0.1, 0.9, 0.19*yf, 0.21*yf, 0.7, 0.9, 0.1, 0.1, 0.0);
      //(void) add_block_of_particles (&pts, 40000, 0.1, 0.9, 0.49*yf, 0.51*yf, 0.1, 0.9, 0.7, 1.0, 0.0);
      //(void) add_block_of_particles (&pts, 10000, 0.1, 0.9, 0.79*yf, 0.81*yf, 0.8, 0.1, 0.8, 10.0, 0.0);

      // set particle image resolution
      if (pnx < 0) pnx = nx;
      if (pny < 0) pny = ny;

      if (use_COLOR) {
         // generate the temporary color image for splatting the particles
         pc[0] = allocate_2d_array_f(pnx,pny);
         pc[1] = allocate_2d_array_f(pnx,pny);
         pc[2] = allocate_2d_array_f(pnx,pny);
         // make sure to zero this
         for (int ix=0; ix<pnx; ix++) {
            for (int iy=0; iy<pny; iy++) {
               pc[0][ix][iy] = 0.0;
               pc[1][ix][iy] = 0.0;
               pc[2][ix][iy] = 0.0;
            }
         }
      }
   }

   // -----------------------------------------------
   // Iterate through time

   step = 0;
   while (true) {

      // fprintf(stdout,"\nBegin step %d\n",step);

      if (use_PARTICLES) {
         // how many particles should there be?
         int npart = 0;
         if (step < part_add_step_start) {
            npart = 0;
         } else if (step < part_add_step_end) {
            npart = (int)(part_add_count_end*
                    (float)(step-part_add_step_start)/(float)(part_add_step_end-part_add_step_start));
         } else if (step < part_rem_step_start) {
            npart = part_add_count_end;
         } else if (step < part_rem_step_end) {
            npart = part_rem_count_end + (int)((part_add_count_end-part_rem_count_end)*
                    (float)(part_rem_step_end-step)/(float)(part_rem_step_end-part_rem_step_start));
         } else {
            npart = part_rem_count_end;
         }

         if (npart > pts.n) {
            // add as many as needed to hit the number
            float thiscol[3];
            const int nadd = npart - pts.n;
            for (int i=0; i<nadd; ++i) {
               if (use_color_linear || use_color_area) {
                  (void) get_color (cmi, cnx, cny, rand()/(float)RAND_MAX, rand()/(float)RAND_MAX, thiscol);
               } else {
                  const float grey = rand()/(float)RAND_MAX;
                  thiscol[0] = grey; thiscol[1] = grey; thiscol[2] = grey;
               }
               const float brite = 0.3*thiscol[0] + 0.6*thiscol[1] + 0.1*thiscol[2];
               (void) add_one_particle (&pts, rand()/(float)RAND_MAX,
                                              yf*rand()/(float)RAND_MAX,
                                              0.0, 0.0,
                                              thiscol[0], thiscol[1], thiscol[2],
                                              expf(1.5*brite)-0.8, part_ballistic);
            }
         } else if (npart < pts.n) {
           // remove some particles
           pts.n = npart;
           if (pts.n < 0) pts.n = 0;
         }
      }

      // ----------------------------
      // write output
      if (writeOutput && step%writeevery == 0) {
      //#pragma omp parallel sections private(outfileroot)
      {
         //#pragma omp section
         if (print_vort) {
            sprintf(outfileroot,"vort_%06d",step);
            write_png (outfileroot,nx,ny,false,use_16bpp,
                       a[W2],-vortscale,2.*vortscale,
                       NULL,0.0,1.0,
                       NULL,0.0,1.0);
         }
         //#pragma omp section
         if (print_mask && use_MASK) {
            sprintf(outfileroot,"mask_%06d",step);
            write_png (outfileroot,nx,ny,false,false,
                       mask,0.0,1.0,
                       NULL,0.0,1.0,
                       NULL,0.0,1.0);
         }
         //#pragma omp section
         if (print_temp) {
            sprintf(outfileroot,"temp_%06d",step);
            if (use_color_linear) {
               //float tempramp = 60.0;
               float mult = 1.0;
               // map temp to colors
               for (int ix=0; ix<nx; ix++) {
                  for (int iy=0; iy<ny; iy++) {
                     float fval = 0.5 * (a[SF][ix][iy] + 1.0);
                     if (fval < 0.0) fval = 0.0;
                     if (fval > 1.0) fval = 1.0;
                     fval = 0.5 + fval*(float)(cmn-1);
                     int ival = floor(fval);
                     float frem = fval - (float)ival;
                     //mult = 1.0 - (1.0-abs(a[SF][ix][iy])) * (1.0-MIN(1.0, simtime/60.0));
                     c[0][ix][iy] = mult * ((1.0-frem)*cml[0][ival] + frem*cml[0][ival+1]);
                     c[1][ix][iy] = mult * ((1.0-frem)*cml[1][ival] + frem*cml[1][ival+1]);
                     c[2][ix][iy] = mult * ((1.0-frem)*cml[2][ival] + frem*cml[2][ival+1]);
                  }
               }
               // and write the image
               write_png (outfileroot,nx,ny,true,use_16bpp,
                          c[0],0.0,1.0,
                          c[1],0.0,1.0,
                          c[2],0.0,1.0);
            } else {
               write_png (outfileroot,nx,ny,false,false,
                          a[SF],-1.0,2.0,
                          NULL,0.0,1.0,
                          NULL,0.0,1.0);
            }
         }
         //#pragma omp section
         if (use_COLOR) {
            if (use_PARTICLES) {
               // convert particle fade frames to a factor here
               const float factor = expf(-1.0/(float)part_draw_fade_frames);
               // merge color field and particle splats
               draw_particles(&pts, yf,
                              nx, ny, 0.0, a[RR], a[GG], a[BB],
                              pnx, pny, factor, pc[0], pc[1], pc[2],
                              part_draw_fac, part_draw_mass_pow, part_draw_vel_pow);
               // draw over it with the mask
               if (use_MASK) {
                  (void) mult_part_by_mask (mask_color_mult, yf, nx,ny,mask, pnx,pny,pc[0],pc[1],pc[2]);
               }
               // and write the png
               sprintf(outfileroot,"out_%06d",step);
               write_png (outfileroot,pnx,pny,true,use_16bpp,
                          pc[0],0.0,1.0,
                          pc[1],0.0,1.0,
                          pc[2],0.0,1.0);
               // now optionally diffuse the particle colors (before we draw again)
               if (part_img_diffus > 0.0) {
                  (void) diffuse_color_img (part_img_diffus, pnx, pny, xbdry, ybdry, pc[0],pc[1],pc[2]);
               }
            } else {
               sprintf(outfileroot,"out_%06d",step);
               write_png (outfileroot,nx,ny,true,use_16bpp,
                          a[RR],0.0,1.0,
                          a[GG],0.0,1.0,
                          a[BB],0.0,1.0);
            }
         }
         //#pragma omp section
         if (print_vel) {
            sprintf(outfileroot,"vel_%06d",step);
            write_png (outfileroot,nx,ny,true,use_16bpp,
                       u[XV],-velscale,2.*velscale,
                       u[YV],-velscale,2.*velscale,
                       a[W2],-vortscale,2.*vortscale);
         }
         //#pragma omp section
         if (print_mu) {
            sprintf(outfileroot,"mu_%06d",step);
            write_png (outfileroot,nx,ny,false,use_16bpp,
                       a[MD],mdlow,mdhigh-mdlow,
                       NULL,0.0,1.0,
                       NULL,0.0,1.0);
         }
      } // OMP sections
      } // if writeOutput && step%writeevery == 0

      // compute wmax and write checkpoint files
      if (checkpointevery > 0 && step > 0 && step%checkpointevery == 0) {
         fprintf(stderr,"  writing checkpoint files for restart...\n");

         // find vort max
         float wmax = 1.e-20;
         for (int i=0; i<nx; i++) {
            for (int j=0; j<ny; j++) {
               if (fabs(a[W2][i][j]) > wmax) wmax = a[W2][i][j];
            }
         }

         // write vorticity (always 16bpp to keep precision)
         strcpy(outfileroot,"restart_vort");
         write_png (outfileroot,nx,ny,false,use_16bpp,
                    a[W2],-wmax,2.*wmax,
                    NULL,0.0,1.0,
                    NULL,0.0,1.0);

         // write temperature
         if (use_TEMP) {
            strcpy(outfileroot,"restart_temp");
            write_png (outfileroot,nx,ny,false,use_16bpp,
                       a[SF],-1.0,2.0,
                       NULL,0.0,1.0,
                       NULL,0.0,1.0);
         }

         // write color
         if (use_COLOR) {
            strcpy(outfileroot,"restart_color");
            write_png (outfileroot,nx,ny,true,use_16bpp,
                       a[RR],0.0,1.0,
                       a[GG],0.0,1.0,
                       a[BB],0.0,1.0);
         }

         // echo restart command line
         fprintf(stderr,"  restart from this point with:\n");
         fprintf(stderr,"%s",argv[0]);
         for (int i=1; i<argc; i++) {
            fprintf(stderr," %s",argv[i]);
         }
         fprintf(stderr," -vf restart_vort.png");
         fprintf(stderr," -vscale %g",wmax);
         if (use_TEMP) fprintf(stderr," -tf restart_temp.png");
         if (use_COLOR) fprintf(stderr," -cf restart_color.png");
         fprintf(stderr,"\n");
      }

      // check end conditions
      if (step >= maxstep) break;

      // accept input (only if interactive scheme)

      // adjust diffusion coefficients - decay past step 300
      if (false && step > 300) {
        for (int i=0; i<MAX_SCALARS; i++) {
          // this will work for the negative diffusivity trigger, too
          sc_diffus[i] *= 1.1;
        }
      }

      // reset freestream
      if (dynamic_freestream) {
         if (step < fs_start_step) {
            freestream[0] = fs_start[0];
            freestream[1] = fs_start[1];
         } else if (step < fs_end_step) {
            const float factor = (float)(step-fs_start_step) / (float)(fs_end_step-fs_start_step);
            freestream[0] = fs_start[0] + factor*(fs_end[0]-fs_start[0]);
            freestream[1] = fs_start[1] + factor*(fs_end[1]-fs_start[1]);
         } else {
            freestream[0] = fs_end[0];
            freestream[1] = fs_end[1];
         }
      }

      // if we want a constant cn time step, compute it here
      if (courantconst > 0.0) {
         // average vmax over a number of steps
         float vlast = ccnvmax[ccnidx];
         ccnidx = (ccnidx+1)%VMAXAVG;
         ccnvmax[ccnidx] = vmax;
         // if any were zero, reset them
         for (int i=0; i<VMAXAVG; i++) if (ccnvmax[i] < 1.e-6) ccnvmax[i] = vmax;
         // find new thing
         float vavg = 0.0;
         for (int i=0; i<VMAXAVG; i++) vavg += ccnvmax[i];
         vavg /= (float)VMAXAVG;
         // recalculate dt
         float newdt = (courantconst/vavg)/(nx+1);
         // and test vs. extrapolation of future vmax
         float testdt = (courantconst/(2.0*vmax-vlast))/(nx+1);
         if (testdt < newdt) newdt = testdt;
         // do not let dt change too much!
         if (newdt > 1.4*dt) dt *= 1.4;
         else if (newdt < 0.5*dt) dt *= 0.5;
         else dt = newdt;
         if (!silent) fprintf(stderr,"  changing dt to %g\n",dt);
      }

      // set the timer
      gettimeofday(&t_last, 0);

      // should we test for freeze?
      if (freeze_flow) {
         if (step > freeze_at) {
            recalc_vel = false;
            move_colors = true;
         } else {
            recalc_vel = true;
            move_colors = false;
         }
      }

      // should we test for stop?
      if (stop_flow) {
         if (step > stop_at) {
            recalc_vel = false;
            move_colors = false;
         } else {
            recalc_vel = true;
            move_colors = true;
         }
      }

      // take one computational step forward in time
      int numsubsteps = 1;
      if (use_MASK) {
         numsubsteps = 1 + (int)(vmax*(nx+1)*dt/10.);
         fprintf(stderr,"need %d substeps\n", numsubsteps); fflush(stderr);
      }
      for (int istep=0; istep<numsubsteps; istep++) {
         effective_re = step_forward_2d (silent,step,isStam,4,
                           nx,ny,xbdry,ybdry,freestream,wallvel,
                           recalc_vel,move_colors,
                           u,a,t,use_MASK,mask,maskerr,sc_diffus,
                           gravtype,gravity,(dt/(float)numsubsteps),
                           use_strong_strat,bn,dens_ratio,acc);
      }

      // move particles
      if (use_PARTICLES) {
         (void) move_particles (&pts, nx,ny,xbdry,ybdry,yf, mask,u[XV],u[YV],a[SF],gravity,particle_speed*dt);
      }

      // read in the image again and overlay it!
      if (use_COLOR && reread_color) {
         // read color PNG of exactly nx by ny resolution into
         //   the given fields and overlay it!
         if (use_color_img && reread_color) {
            //read_png(colorfilename,nx,ny,true,true,0.01,true,
            //   a[RR],0.0,1.0,a[GG],0.0,1.0,a[BB],0.0,1.0);
            read_png(colorfilename,nx,ny,true,
                     true,overlay_fraction,darkenonly,
                     a[RR],0.0,1.0,a[GG],0.0,1.0,a[BB],0.0,1.0);
            // use the brightness component to indicate density!
            //for (int ix=0; ix<nx; ix++) {
            //   for (int iy=0; iy<ny; iy++) {
            //      a[SF][ix][iy] = 1.0 - 0.6*a[GG][ix][iy] -
            //                       0.3*a[RR][ix][iy] -
            //                       0.1*a[BB][ix][iy];
            //   }
            //}
         }
      }

      // read in the vorticity image and overlay it
      if (recalc_vel && vort_reread && use_vort_img) {
         // read grayscale PNG of exactly nx by ny resolution
         // the 2 is to force the new data to be simply added to the old
         // the funny min bounds are to allow value of 127 to become 0.0
         read_png(vortfilename,nx,ny,false,
                  2,1.0,false,
                  a[W2],-254.*vortscale/255.,2.0*vortscale,
                  NULL, -254.*vortscale/255.,2.0*vortscale,
                  NULL, -254.*vortscale/255.,2.0*vortscale);
      }

      if (false) {
         // open or close blocks in the mask
         (void) update_mask_with_blocks_2 (mask, a, nx, ny, simtime, dt);
         // optionally generate repeatedly-overlaid mask
         (void) overlay_mask (nx, ny, mask);
      }

      if (dynamic_mask) {
         (void) set_mask_from_temporal (step, nx, ny, mask, maskfilename,
                                        dmask_step_start, dmask_val_start,
                                        dmask_step_end, dmask_val_end,
                                        dmask_width, dmask_power);
      }

      if (false && step%50 == 0 && step < 1001) {
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
         for (int ix=0; ix<nx; ix++) {
            for (int iy=0; iy<ny; iy++) {
               a[SF][ix][iy] = 1.0 - 0.6*a[GG][ix][iy] -
                                0.3*a[RR][ix][iy] -
                                0.1*a[BB][ix][iy];
            }
         }
      }

      if (false && use_MASK) {
         // allocate space
         if (first_time) {
            shear = allocate_2d_array_f(nx,ny);
         }

         // adjust color value based on shear present in flow
         (void) find_shear_magnitude (nx, ny, xbdry, ybdry, u[XV],1.0, u[YV],1.0, shear);

         float minshear = 9.9e+9;
         float maxshear = -9.9e+9;
         float meanshear = 0.0;

         //#pragma omp parallel for private() reduce()
         for (int ix=0; ix<nx; ix++) {
            for (int iy=0; iy<ny; iy++) {
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

         if (true) {
            sprintf(outfileroot,"shear_%06d",step);
            write_png (outfileroot,nx,ny,false,false,
                       shear,-100.0,200.0,
                       NULL,0.0,1.0,
                       NULL,0.0,1.0);
         }
      }

      // use vorticity and mask to modify mask
      if (false && use_MASK) {
         (void) mod_mask_with_vort(nx,ny,dt,mask,a[W2]);
      }

      // use velocity and mask to modify mask
      if (false && use_MASK) {
         (void) mod_mask_with_vel(nx,ny,dt,mask,u[XV],u[YV]);
      }

      // overlay more random vorticity
      if (vort_add_rand) {
         // add a random field of vorticity over the existing vorticity
         float factor = 2.0 * vort_add_factor * dt;
         bool doit = false;
         if (step < vort_add_start) {
            doit = false;
         } else if (step < vort_add_peak) {
            doit = true;
            // lerp magnitude
            factor *= (float)(step-vort_add_start)/(float)(vort_add_peak-vort_add_start);
         } else if (step < vort_add_down) {
            doit = true;
            // constant peak magnitude
         } else if (step < vort_add_end) {
            doit = true;
            // lerp magnitude
            factor *= (float)(vort_add_end-step)/(float)(vort_add_end-vort_add_down);
         } else {
            doit = false;
         }
         if (doit) {
            #pragma omp parallel for
            for (int ix=0; ix<nx; ix++) {
               for (int iy=0; iy<ny; iy++) {
                  a[W2][ix][iy] += factor * (rand()/(float)RAND_MAX - 0.5);
               }
            }
         }
      }

      // supppress vorticity according to the square
      if (vort_suppress) {
         // multiply entire field by a decay factor
         const float factor = vort_suppress_factor * dt;
         #pragma omp parallel for
         for (int ix=0; ix<nx; ix++) {
            for (int iy=0; iy<ny; iy++) {
               a[W2][ix][iy] *= exp(-factor*a[W2][ix][iy]*a[W2][ix][iy]);
            }
         }
      }

      // clamp the temperature field
      if (use_TEMP && true) {
         #pragma omp parallel for
         for (int ix=0; ix<nx; ix++) {
            for (int iy=0; iy<ny; iy++) {
               if (a[SF][ix][iy] > 1.0) a[SF][ix][iy] = 1.0;
               if (a[SF][ix][iy] < -1.0) a[SF][ix][iy] = -1.0;
            }
         }
      }

      // overlay the heat map
      if (use_heat_img) {
         #pragma omp parallel for
         for (int ix=0; ix<nx; ix++) {
            for (int iy=0; iy<ny; iy++) {
               a[SF][ix][iy] += dt * heat_coeff * heat[ix][iy];
               if (a[SF][ix][iy] > 1.0) a[SF][ix][iy] = 1.0;
               if (a[SF][ix][iy] < -1.0) a[SF][ix][iy] = -1.0;
            }
         }
      }

      // overlay the edge color input
      if (use_COLOR && overlay_color_tb) {
         int nthick = 6;
         for (int iy=0; iy<nthick; iy++) {
            float tfac = 1.0 / (float)(iy+1);
            for (int ix=0; ix<nx; ix++) {
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
         for (int ix=0; ix<nthick; ix++) {
            float tfac = 1.0 / (float)(ix+1);
            for (int iy=0; iy<ny; iy++) {
               a[RR][ix][iy] = tfac*color_left[0][iy] + (1.0-tfac)*a[RR][ix][iy];
               a[GG][ix][iy] = tfac*color_left[1][iy] + (1.0-tfac)*a[GG][ix][iy];
               a[BB][ix][iy] = tfac*color_left[2][iy] + (1.0-tfac)*a[BB][ix][iy];
               a[RR][nx-1-ix][iy] = tfac*color_right[0][iy] + (1.0-tfac)*a[RR][nx-1-ix][iy];
               a[GG][nx-1-ix][iy] = tfac*color_right[1][iy] + (1.0-tfac)*a[GG][nx-1-ix][iy];
               a[BB][nx-1-ix][iy] = tfac*color_right[2][iy] + (1.0-tfac)*a[BB][nx-1-ix][iy];
            }
         }
      }

      // nudge colors toward a given color
      if (use_COLOR && use_color_area && false) {
         const float tfac = 0.1;
         float thiscol[3];
         for (int ix=0; ix<nx; ix++) {
            for (int iy=0; iy<ny; iy++) {
               // use current vel to select color
               const float thisvel = sqrt(pow(u[XV][ix][iy],2) + pow(u[YV][ix][iy],2));
               (void) get_color (cmi, cnx, cny, 0.02*simtime, MIN(1.0,5.0*thisvel), thiscol);
               a[RR][ix][iy] += tfac * (thiscol[0] - a[RR][ix][iy]);
               a[GG][ix][iy] += tfac * (thiscol[1] - a[GG][ix][iy]);
               a[BB][ix][iy] += tfac * (thiscol[2] - a[BB][ix][iy]);
            }
         }
      }
      if (use_COLOR && use_cint) {
         if (step%cint_interval == cint_interval-1) {
            // scale the colors
            float factor = 1.0f / (float)cint_interval;
            for (int ix=0; ix<nx; ix++) { for (int iy=0; iy<ny; iy++) { cint[0][ix][iy] *= factor; } }
            for (int ix=0; ix<nx; ix++) { for (int iy=0; iy<ny; iy++) { cint[1][ix][iy] *= factor; } }
            for (int ix=0; ix<nx; ix++) { for (int iy=0; iy<ny; iy++) { cint[2][ix][iy] *= factor; } }
            // write the colors
            sprintf(cintfileroot,"cint_%06d",step+1);
            write_png (cintfileroot,nx,ny,true,use_16bpp,
                       cint[0],0.0,1.0,
                       cint[1],0.0,1.0,
                       cint[2],0.0,1.0);
            // and zero them again
            for (int ix=0; ix<nx; ix++) { for (int iy=0; iy<ny; iy++) { cint[0][ix][iy] = 0.0f; } }
            for (int ix=0; ix<nx; ix++) { for (int iy=0; iy<ny; iy++) { cint[1][ix][iy] = 0.0f; } }
            for (int ix=0; ix<nx; ix++) { for (int iy=0; iy<ny; iy++) { cint[2][ix][iy] = 0.0f; } }
         } else {
            // accumulate the color
            for (int ix=0; ix<nx; ix++) { for (int iy=0; iy<ny; iy++) { cint[0][ix][iy] += a[RR][ix][iy]; } }
            for (int ix=0; ix<nx; ix++) { for (int iy=0; iy<ny; iy++) { cint[1][ix][iy] += a[GG][ix][iy]; } }
            for (int ix=0; ix<nx; ix++) { for (int iy=0; iy<ny; iy++) { cint[2][ix][iy] += a[BB][ix][iy]; } }
         }
      }

      // calculate the total time elapsed, and the time for this last step
      gettimeofday(&t_curr, 0);
      walltime = (t_curr.tv_sec - t_last.tv_sec) + (t_curr.tv_usec - t_last.tv_usec)*1e-6;

      // write statistics
      if (writeOutput) {
         vmax = compute_and_write_stats(silent,step,dt,simtime,walltime,nx,ny,xbdry,ybdry,u,a[W2]);
      }

      if (first_time) first_time = false;

      // ----------------------------
      step++;
      simtime += dt;
   }

   // close out the stats file properly
   if (writeOutput && maxstep > 0) {
      compute_and_write_stats(silent,-1,dt,simtime,walltime,nx,ny,xbdry,ybdry,u,a[W2]);
   }

   if (!silent) fprintf(stderr,"\nDone.\n");
   exit(0);
}


/*
 * Do just what it says
 */
float compute_and_write_stats(int silent, int step, float dt, float simtime, double walltime,
      int nx, int ny, int xbdry, int ybdry, float ***u, float **w) {

   float ke,base_mult,multx,multy,vmax,wmax,cn;
   static double total_time = 0.0;
   static char outfile[MAXCHARS];
   static bool initialized = false;
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
      initialized = true;
      fprintf(outp,"# step, sim time, KE, vort max, vel max, CN, cpu time step, cpu time total\n");
      fflush(outp);
   }

   // find kinetic energy
   ke = 0.;
   base_mult = 1.0 / ((nx-1) * (ny-1));
   for (int i=0; i<nx; i++) {
      if (i==0 || i==nx-1) multx = base_mult*0.5;
         else multx = base_mult;
      for (int j=0; j<ny; j++) {
         if (j==0 || j==ny-1) multy = multx*0.5;
            else multy = multx;
         ke += multy*(u[XV][i][j]*u[XV][i][j] + u[YV][i][j]*u[YV][i][j]);
      }
   }
   ke /= 2.;
   //fprintf(stdout,"ke is %g\n",ke);

   // find vort max
   wmax = 0.;
   for (int i=0; i<nx; i++) {
      for (int j=0; j<ny; j++) {
         if (fabs(w[i][j]) > wmax) wmax = w[i][j];
      }
   }

   // find vmax
   vmax = 0.;
   for (int i=0; i<nx; i++) {
      for (int j=0; j<ny; j++) {
         float velsq = u[XV][i][j]*u[XV][i][j] + u[YV][i][j]*u[YV][i][j];
         if (velsq > vmax) vmax = velsq;
      }
   }
   vmax = sqrt(vmax);
   cn = vmax*(nx+1)*dt;

   // and time
   total_time += walltime;

   // write the stats here
   if (!silent) {
      if (step%10 == 0) fprintf(stdout,"# step, sim time, KE, vort max, vel max, CN, cpu step time, cpu total time\n");
      fprintf(stdout,"%d %g %g %g %g %g %g %g\n",step,simtime,ke,wmax,vmax,cn,walltime,total_time);
      fflush(stdout);
   }
   fprintf(outp,"%d %g %g %g %g %g %g %g\n",step,simtime,ke,wmax,vmax,cn,walltime,total_time);
   fflush(outp);

   return(vmax);
}


/*
 * Drop a groovy paint splat
 */
void paint_splat (float sx,float sy, float ex,float ey, float rad,float thresh,
                  float r,float g,float b, int nspots,
                  int nx,int ny, float **red,float **grn,float **blu) {

   float cx,cy,px,py,temp;
   float **spot;

   // must use isosurface of Gaussian-splattered spots
   // first, create x,y,z locations and strengths of blobs
   spot = allocate_2d_array_f(nspots,3);
   for (int i=0; i<nspots; i++) {
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
   for (int ix=0; ix<nx; ix++) {
      for (int iy=0; iy<ny; iy++) {
         // pixel position
         px = ix/(float)nx;
         py = iy/(float)ny;
         // scalar potential
         temp = 0.;
         for (int i=0; i<nspots; i++) {
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
void get_random_color (float ***cmi, int nx, int ny, float *thisc) {

   // choose random coordinate
   int ix = nx*(float)rand()/(float)RAND_MAX;
   int iy = ny*(float)rand()/(float)RAND_MAX;

   thisc[0] = cmi[0][ix][iy];
   thisc[1] = cmi[1][ix][iy];
   thisc[2] = cmi[2][ix][iy];
}

/*
 * grab a color from the color image
 */
void get_color (float ***cmi, int imx, int imy, float x, float y, float *thisc) {

   const int ix = x * imx;
   const int iy = y * imy;

   if (ix>-1 && iy>-1 && ix<imx && iy<imy) {
     thisc[0] = cmi[0][ix][iy];
     thisc[1] = cmi[1][ix][iy];
     thisc[2] = cmi[2][ix][iy];
   } else {
     thisc[0] = cmi[0][0][0];
     thisc[1] = cmi[1][0][0];
     thisc[2] = cmi[2][0][0];
   }
}


/*
 * This function writes basic usage information to stderr,
 * and then quits. Too bad.
 */
int Usage(char progname[MAXCHARS],int status) {

   static char **cpp, *help_message[] = {
   "where [-options] are one or more of the following:                         ",
   "                                                                           ",
   "   -x [int]    resolution in x-direction; default=512                      ",
   "   -y [int]    resolution in y-direction; default=512                      ",
   "                                                                           ",
   "   -px         make domain periodic in x-direction; default is wall bc     ",
   "   -py         make domain periodic in y-direction; default is wall bc     ",
   "   -open       make domain open in both directions; default is wall bc     ",
   "   -bcl [float] set wall-tangent boundary vel on left side                  ",
   "   -bcr [float] set wall-tangent boundary vel on right side                 ",
   "   -bcb [float] set wall-tangent boundary vel on bottom side                ",
   "   -bct [float] set wall-tangent boundary vel on top side                   ",
   "                                                                           ",
   "   -fs [float] [float] set freestream in x and y directions                ",
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
   "   -pa [stepstart] [stepend] [endcount]                                    ",
   "               add particles, at stepstart there should be no particles,   ",
   "               and at stepend there should be endcount                     ",
   "   -pr [stepstart] [stepend] [endcount]                                    ",
   "               remove particles starting at step stepstart, at step        ",
   "               stepend there should be endcount                            ",
   "   -ps [float] multiplier on particle speed                                ",
   "   -pb [float] multiplier on particle ballistic coefficient                ",
   "   -pf [float] number of frames to fade particles by factor of 1/e         ",
   "   -pd [factor] [mass power] [vel power]                                   ",
   "               scale particle drawing value using these properties         ",
   "                                                                           ",
   "   -va [float] add random vorticity every time step, value is magnitude    ",
   "   -vr [float] suppress vorticity according to its square, val is coeff    ",
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
   "   -checkpoint [int] write checkpoint files every so many steps,           ",
   "               default=-1 (no checkpointing)                               ",
   "   -cint [int] integrate colors over a number of steps, default=off        ",
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

