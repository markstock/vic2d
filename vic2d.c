/*
 * vic2d
 *
 * Copyright 2004-8 Mark J. Stock mstock@umich.edu
 *
 * a 2D vortex method which uses the method of characteristics for the
 * convection step, and a single explicit step for diffusion and vorticity
 * creation from walls
 *
 */

#include "vicmoc.h"
void compute_and_write_stats(int, float, float, float, int, int, int, int, float***);

int main(int argc,char **argv) {

   int isStam = FALSE;			// set to TRUE to run Stable Fluids calc
   int i,j;
   int nx,ny;
   int xbdry,ybdry;
   int ix,iy;
   int step = -1;
   int maxstep = 0;
   int outscale;

   int sc_cnt = 0;
   int sc_type[MAX_SCALARS];
   float sc_diffus[MAX_SCALARS];
   int use_TEMP = FALSE;
   int use_DF = FALSE;
   int use_COLOR = FALSE;
   int use_MASK = FALSE;
   int use_image = FALSE;
   int sfi,red,green,blue;

   int writeOutput = TRUE;
   int print_vort = FALSE;
   int print_vel = FALSE;

   float courant;			// non-dim time step size
   float reynolds;			// non-dim viscosity
   float prandtl;			// non-dim thermal diffusivity
   float schmidt;			// non-dim scalar diffusivity
   float bn;				// Boussinesq coefficient
   float dt;
   float *freestream;
   int gravtype;
   float *gravity;
   float coeff,yf,thickness;
   float effective_re;
   float px,py,cx,cy,rad,temp,temp2,maxval,minval,scale;
   float vortscale,velscale;
   float cputime = 0.;
   float totcputime = 0.;
   float avgcputime = 0.;
   float simtime = 0.;

   float **u[2];			// velocities
   float **a[MAX_SCALARS];		// other scalars
   float **t[MAX_SCALARS];		// temporary velocities, vorticity, etc
   float **mask;			// flow mask
   float **heat;			// constant heat source map

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


   // set default simulation properties
   nx = 129;
   ny = 129;
   xbdry = WALL;
   ybdry = WALL;
   courant = 10.0;
   dt = -1.0;
   reynolds = 1000.0;
   prandtl = 0.7;
   schmidt = 1.;
   outscale = 1;
   maxstep = 10000;
   bn = 1.;
   for (i=0; i<MAX_SCALARS; i++) sc_type[i] = -1;
   writeOutput = TRUE;
   freestream = allocate_1d_array_f(2);
   freestream[0] = 0.0;
   freestream[1] = 0.0;
   gravtype = 0;
   gravity = allocate_1d_array_f(2);
   gravity[0] = 0.0;
   gravity[1] = 1.0;
   vortscale = 10.;
   velscale = 1.;

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
      } else if (strncmp(argv[i], "-cn", 3) == 0) {
         courant = atof(argv[++i]);

      } else if (strncmp(argv[i], "-tf", 3) == 0) {
         strcpy (tempfilename,argv[++i]);
         use_temp_img = TRUE;
         use_TEMP = TRUE;
      } else if (strncmp(argv[i], "-heat", 4) == 0) {
         strcpy (heatfilename,argv[++i]);
         use_heat_img = TRUE;
         use_TEMP = TRUE;
      } else if (strncmp(argv[i], "-t", 2) == 0) {
         use_TEMP = TRUE;
      } else if (strncmp(argv[i], "-b", 2) == 0) {
         bn = atof(argv[++i]);
      } else if (strncmp(argv[i], "-pr", 3) == 0) {
         prandtl = atof(argv[++i]);

      } else if (strncmp(argv[i], "-cf", 3) == 0) {
         strcpy (colorfilename,argv[++i]);
         use_color_img = TRUE;
         use_COLOR = TRUE;
      } else if (strncmp(argv[i], "-c", 2) == 0) {
         use_COLOR = TRUE;
      } else if (strncmp(argv[i], "-sc", 3) == 0) {
         schmidt = atof(argv[++i]);

      } else if (strncmp(argv[i], "-vf", 3) == 0) {
         strcpy (vortfilename,argv[++i]);
         use_vort_img = TRUE;
      } else if (strncmp(argv[i], "-re", 3) == 0) {
         reynolds = atof(argv[++i]);

      } else if (strncmp(argv[i], "-df", 3) == 0) {
         strcpy (divfilename,argv[++i]);
         use_div_img = TRUE;
      } else if (strncmp(argv[i], "-div", 3) == 0) {
         use_DF = TRUE;

      } else if (strncmp(argv[i], "-mf", 3) == 0) {
         strcpy (maskfilename,argv[++i]);
         use_mask_img = TRUE;
         use_MASK = TRUE;
      } else if (strncmp(argv[i], "-m", 2) == 0) {
         use_MASK = TRUE;

      } else if (strncmp(argv[i], "-os", 3) == 0) {
         outscale = atoi(argv[++i]);
      } else if (strncmp(argv[i], "-dt", 3) == 0) {
         dt = atof(argv[++i]);
      } else if (strncmp(argv[i], "-step", 5) == 0) {
         maxstep = atoi(argv[++i]);
      } else if (strncmp(argv[i], "-stam", 5) == 0) {
         isStam = TRUE;
      } else if (strncmp(argv[i], "-fs", 3) == 0) {
         freestream[0] = atof(argv[++i]);
         freestream[1] = atof(argv[++i]);
         xbdry = OPEN;
         ybdry = OPEN;

      } else if (strncmp(argv[i], "-pvort", 4) == 0) {
         print_vort = TRUE;
      } else if (strncmp(argv[i], "-vortscale", 6) == 0) {
         vortscale = atof(argv[++i]);
      } else if (strncmp(argv[i], "-pvel", 4) == 0) {
         print_vel = TRUE;
      } else if (strncmp(argv[i], "-velscale", 5) == 0) {
         velscale = atof(argv[++i]);
      } else if (strncmp(argv[i], "-noout", 6) == 0) {
         writeOutput = FALSE;
      } else {
         fprintf(stderr,"Unknown option (%s)\n",argv[i]);
         (void) Usage(progname,0);
      }
   }

   // set y max bound (if square, yf=1.0, otherwise scale it as xmax=1.0
   yf = (float)(ny-1)/(float)(nx-1);

   // set time step (dt), base it on the maximum courant number (max speed of)
   //  but only if courant number is not set
   if (dt < 0.)
     dt = courant/(nx-1);
   else
     courant = dt*(nx-1);
   fprintf(stdout,"Using dt=%g, or courant number=%g\n",dt,courant);

   // -----------------------------------------------
   // Allocate and initialize

   // velocities first
   u[XV] = allocate_2d_array_f(nx,ny);
   u[YV] = allocate_2d_array_f(nx,ny);

   // initialize the scalar count
   sc_cnt = 0;

   // one necessary scalar type is vorticity, always at position 0
   a[sc_cnt] = allocate_2d_array_f(nx,ny);
   t[sc_cnt] = allocate_2d_array_f(nx,ny);
   sc_type[sc_cnt] = W2;
   // set coefficient of momentum diffusivity
   sc_diffus[sc_cnt] = 1.0/reynolds;
   sc_cnt++;

   // split on solution type
   if (isStam) {

      // both velocities are scalars, as they must be diffused normally
      a[sc_cnt] = allocate_2d_array_f(nx,ny);
      t[sc_cnt] = allocate_2d_array_f(nx,ny);
      sc_type[sc_cnt] = XV;
      // set coefficient of momentum diffusivity
      sc_diffus[sc_cnt] = 1.0/reynolds;
      sc_cnt++;

      a[sc_cnt] = allocate_2d_array_f(nx,ny);
      t[sc_cnt] = allocate_2d_array_f(nx,ny);
      sc_type[sc_cnt] = YV;
      sc_diffus[sc_cnt] = 1.0/reynolds;
      sc_cnt++;
   }

   // set more scalars; this one: density (or temperature)
   if (use_TEMP) {
      sc_type[sc_cnt] = SF;
      // set thermal diffusivity
      sc_diffus[sc_cnt] = 1./(reynolds*prandtl);
      a[sc_cnt] = allocate_2d_array_f(nx,ny);
      t[sc_cnt] = allocate_2d_array_f(nx,ny);
      sc_cnt++;
   }

   // this one: divergence
   if (use_DF) {
      sc_type[sc_cnt] = DF;
      // set dilation diffusivity equal to momentum diffusivity,
      // but should really be second viscosity (which is equivalent or greater
      // than momentum diffusivity mu)
      sc_diffus[sc_cnt] = 1.0/reynolds;
      a[sc_cnt] = allocate_2d_array_f(nx,ny);
      t[sc_cnt] = allocate_2d_array_f(nx,ny);
      sc_cnt++;
   }

   // this one: color
   if (use_COLOR) {
      sc_type[sc_cnt] = RR;
      // set scalar diffusivity
      sc_diffus[sc_cnt] = 1./(reynolds*schmidt);
      a[sc_cnt] = allocate_2d_array_f(nx,ny);
      t[sc_cnt] = allocate_2d_array_f(nx,ny);
      sc_cnt++;
      // now green
      sc_type[sc_cnt] = GG;
      sc_diffus[sc_cnt] = 1./(reynolds*schmidt);
      a[sc_cnt] = allocate_2d_array_f(nx,ny);
      t[sc_cnt] = allocate_2d_array_f(nx,ny);
      sc_cnt++;
      // now blue
      sc_type[sc_cnt] = BB;
      sc_diffus[sc_cnt] = 1./(reynolds*schmidt);
      a[sc_cnt] = allocate_2d_array_f(nx,ny);
      t[sc_cnt] = allocate_2d_array_f(nx,ny);
      sc_cnt++;
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
      read_png(vortfilename,nx,ny,FALSE,FALSE,1.0,
         a[W2],-1.0,2.0,NULL,-1.0,2.0,NULL,-1.0,2.0);
   }


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
      if (TRUE) {
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
         read_png(maskfilename,nx,ny,FALSE,FALSE,1.0,
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
      for (ix=0; ix<nx; ix++) {
         for (iy=0; iy<ny; iy++) {
            mask[ix][iy] = (maxval-mask[ix][iy])/temp;
            if (mask[ix][iy] < 0.0) mask[ix][iy] = 0.0;
            if (mask[ix][iy] > 1.0) mask[ix][iy] = 1.0;
         }
      }
   }

   // Initial velocity solve -------------------------

   // if you want to show velocity on the first step, solve for it here
   find_vels_2d (step,isStam,nx,ny,xbdry,ybdry,freestream,u[XV],u[YV],a[W2],use_MASK,mask);


   // set color ---------------------------------------

   // first, find what index the scalar uses
   if (use_COLOR) {
      for (i=0; i<sc_cnt; i++) if (sc_type[i] == RR) red = i;
      for (i=0; i<sc_cnt; i++) if (sc_type[i] == GG) green = i;
      for (i=0; i<sc_cnt; i++) if (sc_type[i] == BB) blue = i;
      // zero the arrays
      for (ix=0; ix<nx; ix++) for (iy=0; iy<ny; iy++) a[red][ix][iy] = 0.0;
      for (ix=0; ix<nx; ix++) for (iy=0; iy<ny; iy++) a[green][ix][iy] = 0.0;
      for (ix=0; ix<nx; ix++) for (iy=0; iy<ny; iy++) a[blue][ix][iy] = 0.0;
      if (use_color_img) {
         // read color PNG of exactly nx by ny resolution into
         //   the given fields
         read_png(colorfilename,nx,ny,TRUE,FALSE,1.0,
            a[red],0.0,1.0,a[green],0.0,1.0,a[blue],0.0,1.0);
      }
      if (FALSE) {
         // left is red, right is green
         for (ix=0; ix<nx/2; ix++) {
            for (iy=0; iy<ny; iy++) {
               a[red][ix][iy] = 1.0;
               a[green][ix][iy] = 0.91372549;
               a[blue][ix][iy] = 0.168627451;
            }
         }
         for (ix=nx/2; ix<nx; ix++) {
            for (iy=0; iy<ny; iy++) {
               a[red][ix][iy] = 0.1171875;
               a[green][ix][iy] = 0.0546875;
               a[blue][ix][iy] = 0.7265625;
            }
         }
      }
   }

   // Set scalar/temperature ----------------------------------

   // first, find what index the scalar uses
   if (use_TEMP) {
      for (i=0; i<sc_cnt; i++) if (sc_type[i] == SF) sfi = i;

      // zero the array
      for (ix=0; ix<nx; ix++) {
         for (iy=0; iy<ny; iy++) {
            a[sfi][ix][iy] = 0.0;
         }
      }

      if (use_temp_img) {
         // read grayscale PNG of exactly nx by ny resolution
         read_png(tempfilename,nx,ny,FALSE,TRUE,1.0,
            a[sfi],0.0,1.0,NULL,0.0,1.0,NULL,0.0,1.0);
      }

      if (use_heat_img) {
         // read grayscale PNG of exactly nx by ny resolution
         read_png(heatfilename,nx,ny,FALSE,FALSE,1.0,
            heat,0.0,1.0,NULL,0.0,1.0,NULL,0.0,1.0);

         for (ix=0; ix<nx; ix++) {
            for (iy=0; iy<ny; iy++) {
               a[sfi][ix][iy] += heat[ix][iy];
               if (a[sfi][ix][iy] > 1.0) a[sfi][ix][iy] = 1.0;
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
                  a[sfi][ix][iy] = 10.0;
            }
         }
      }
      if (FALSE) {
         // create a smooth blob of scalar
         // add_smooth_circular_blob(nx,ny,xbdry,ybdry,a[sfi],0.7,0.15,0.1,1.0);
         add_smooth_circular_blob(nx,ny,xbdry,ybdry,a[sfi],0.4,0.0,0.07,-1.0);
      }
      if (FALSE) {
         // create bands of scalar
         for (ix=0; ix<nx; ix++) {
            for (iy=0; iy<ny; iy++) {
               //a[sfi][ix][iy] = sin((double)(ix*0.2)) +
               //                 0.2*rand()/RAND_MAX - 0.1;
               //a[sfi][ix][iy] = 5.*sin((double)(ix*0.2));
               // run109/06
               //a[sfi][ix][iy] = 5.*sin(410.*(double)(ix)/nx);
               // run109/07
               a[sfi][ix][iy] = 0.5 + 0.5*sin(180.*(double)(ix)/nx);
            }
         }
      }
      if (FALSE) {
         // create square bands of scalar
         for (ix=0; ix<nx; ix++) {
            for (iy=0; iy<ny; iy++) {
               //a[sfi][ix][iy] = sin((double)(ix*0.2)) +
               //                 0.2*rand()/RAND_MAX - 0.1;
               px = fabs(ix - (float)nx/2.);
               py = fabs(iy - (float)ny/2.);
               if (px > py) temp = px;
               else temp = py;
               if (temp > nx/3.1) temp = -M_PI/2.;
               else temp *= 0.3;
               a[sfi][ix][iy] = 5.*sin(temp);
            }
         }
      }
      if (FALSE) {
         // create a RTI with a sinusoidal interface
         //thickness = 0.001;
         thickness = 2.0/(float)nx;
         for (ix=0; ix<nx; ix++) {
            px = (ix - (float)nx/2.)/(float)nx;
            yf = 0.02*sin(px*9.*M_PI);
            for (iy=0; iy<ny; iy++) {
               //a[sfi][ix][iy] = sin((double)(ix*0.2)) +
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
               a[sfi][ix][iy] = temp;
            }
         }
      }
      if (FALSE) {
        // create a random field of scalar
        scale = 10.;
        for (ix=0; ix<nx; ix++) {
          for (iy=0; iy<ny; iy++) {
            a[sfi][ix][iy] += scale*rand()/RAND_MAX - 0.5*scale + 0.5;
          }
        }
      }
      if (TRUE && use_COLOR) {
         // use the brightness component to indicate density!
         for (ix=0; ix<nx; ix++) {
            for (iy=0; iy<ny; iy++) {
               a[sfi][ix][iy] = 1.0 - 0.6*a[green][ix][iy] -
                                0.3*a[red][ix][iy] -
                                0.1*a[blue][ix][iy];
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
      if (writeOutput) {
      if (print_vort) {
         sprintf(outfileroot,"vort_%05d",step);
         //write_png (outfileroot,nx,ny,FALSE,FALSE,a[W2],-1.0,2.0,
         //write_png (outfileroot,nx,ny,FALSE,FALSE,a[W2],-10.0,20.0,
         //write_png (outfileroot,nx,ny,FALSE,FALSE,a[W2],-250.0,500.0,
         //write_png (outfileroot,nx,ny,FALSE,FALSE,a[W2],-100.0,200.0,
         //write_png (outfileroot,nx,ny,FALSE,FALSE,a[W2],-10.0,20.0,
         write_png (outfileroot,nx,ny,FALSE,FALSE,a[W2],-vortscale,2.*vortscale,
                          NULL,0.0,1.0,NULL,0.0,1.0);
      }
      if (use_TEMP) {
         sprintf(outfileroot,"temp_%05d",step);
         write_png (outfileroot,nx,ny,FALSE,FALSE,a[sfi],0.0,1.0,
                          NULL,0.0,1.0,NULL,0.0,1.0);
      }
      if (use_COLOR) {
         sprintf(outfileroot,"out_%05d",step);
         write_png (outfileroot,nx,ny,TRUE,FALSE,a[red],0.0,1.0,
                          a[green],0.0,1.0,a[blue],0.0,1.0);
      }
      if (print_vel) {
         sprintf(outfileroot,"vel_%05d",step);
         write_png (outfileroot,nx,ny,TRUE,FALSE,
                    u[XV],-velscale,2.*velscale,
                    u[YV],-velscale,2.*velscale,
                    a[W2],-vortscale,2.*vortscale);
      }
      }

      // check end conditions
      if (step >= maxstep) break;

      // accept input (only if interactive scheme)


      // adjust diffusion coefficients - decay past step 300
      if (FALSE && step > 300) {
        for (i=0; i<sc_cnt; i++) {
          sc_diffus[i] *= 1.1;
        }
      }

      // set the timer
      last_tics = clock();

      // take one computational step forward in time
      effective_re = step_forward_2d (step,isStam,nx,ny,xbdry,ybdry,freestream,
                        u,a,t,use_MASK,mask,sc_cnt,sc_type,sc_diffus,
                        gravtype,gravity,dt,bn);

      // read in the image again and overlay it!
      if (use_COLOR && FALSE) {
         // read color PNG of exactly nx by ny resolution into
         //   the given fields and overlay it!
         if (use_image) {
            read_png(colorfilename,nx,ny,TRUE,TRUE,0.1,
               a[red],0.0,1.0,a[green],0.0,1.0,a[blue],0.0,1.0);
            // use the brightness component to indicate density!
            for (ix=0; ix<nx; ix++) {
               for (iy=0; iy<ny; iy++) {
                  a[sfi][ix][iy] = 1.0 - 0.6*a[green][ix][iy] -
                                   0.3*a[red][ix][iy] -
                                   0.1*a[blue][ix][iy];
               }
            }
         }
      }

      // overlay the heat map
      if (use_heat_img) {
         for (ix=0; ix<nx; ix++) {
            for (iy=0; iy<ny; iy++) {
               a[sfi][ix][iy] += heat[ix][iy];
               if (a[sfi][ix][iy] > 1.0) a[sfi][ix][iy] = 1.0;
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
         compute_and_write_stats(step,dt,simtime,cputime,nx,ny,xbdry,ybdry,u);
      }

      // ----------------------------
      step++;
      simtime += dt;
   }

   // close out the stats file properly
   if (writeOutput) {
      compute_and_write_stats(-1,dt,simtime,cputime,nx,ny,xbdry,ybdry,u);
   }

   fprintf(stderr,"\nDone.\n");
   exit(0);
}


/*
 * Do just what it says
 */
void compute_and_write_stats(int step, float dt, float simtime, float cputime,
      int nx, int ny, int xbdry, int ybdry, float ***u) {

   int i,j;
   float ke,base_mult,multx,multy,velsq,vmax,cn;
   static char outfile[80];
   static int initialized = FALSE;
   static FILE *outp;

   // if last time through, close file pointer and RETURN
   if (step == -1) {
      fclose(outp);
      return;
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
   fprintf(stdout,"%d %g %g %g %g %g\n",step,simtime,ke,vmax,cn,cputime);
   fprintf(outp,"%d %g %g %g %g %g\n",step,simtime,ke,vmax,cn,cputime);
   fflush(outp);

   return;
}


/*
 * This function writes basic usage information to stderr,
 * and then quits. Too bad.
 */
int Usage(char progname[80],int status) {

   static char **cpp, *help_message[] = {
   "where [-options] are one or more of the following:                         ",
   "                                                                           ",
   "   -x [int]    resolution in x-direction; default=128                      ",
   "                                                                           ",
   "   -y [int]    resolution in y-direction; default=128                      ",
   "                                                                           ",
   "   -px         make domain periodic in x-direction; default is not periodic",
   "                                                                           ",
   "   -py         make domain periodic in y-direction; default is not periodic",
   "                                                                           ",
   "   -o [int]    set output resolution to num times simulation resolution;   ",
   "               default=1                                                   ",
   "                                                                           ",
   "   -pvel       write velocity field (u_00000.png and v_00000.png)          ",
   "                                                                           ",
   "   -velscale [float]                                                       ",
   "               scale factor for velocity output                            ",
   "                                                                           ",
   "                                                                           ",
   "   -pvort      write vorticity field (wz_00000.png)                        ",
   "                                                                           ",
   "   -vortscale [float]                                                      ",
   "               scale factor for vorticity output                           ",
   "                                                                           ",
   "   -vf [name]  read PNG file and use as initial vorticity field            ",
   "                                                                           ",
   "   -re [float] set Reynolds number, inverse of momentum diffusivity;       ",
   "               default=1000                                                ",
   "                                                                           ",
   "                                                                           ",
   "   -t          track and print temperature/scalar field (dens_00000.png)   ",
   "                                                                           ",
   "   -tf [name]  read PNG file and use as temperature/scalar field           ",
   "                                                                           ",
   "   -b [float]  set the Boussinesq number; default=1                        ",
   "                                                                           ",
   "   -pr [float] set Prandtl number, inverse of thermal diffusivity;         ",
   "               default=0.7                                                 ",
   "                                                                           ",
   "                                                                           ",
   "   -c          track and print color scalars (out_00000.png)               ",
   "                                                                           ",
   "   -cf [name]  read PNG file and use as color scalar field                 ",
   "                                                                           ",
   "   -sc [float] set Schmidt number, inverse of scalar diffusivity;          ",
   "               default=10.0                                                ",
   "                                                                           ",
   "                                                                           ",
   "   -div        track and print dilation field (div_00000.png)              ",
   "                                                                           ",
   "   -df [name]  read PNG file and use as initial dilation/divergence field  ",
   "                                                                           ",
   "                                                                           ",
   "   -mf [name]  read PNG file and use as flow mask, black areas are solid   ",
   "                                                                           ",
   "                                                                           ",
   "   -dt [float] set the time step increment, if unused, Courant number will ",
   "               determine the time step; no default                         ",
   "                                                                           ",
   "   -cn [float] set the Courant number, based on a speed of 1.0; try        ",
   "               10 to 50 for a small simulation (32 to 128) or 50 to 500    ",
   "               for a larger one (256 to 2048); default=10                  ",
   "                                                                           ",
   "   -step [int] run only num steps; default=1000                            ",
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

