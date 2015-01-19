/*
 * VIC-MOC - libvicmoc.c - the library for the 2D routines
 *
 * Copyright 2004-10 Mark J. Stock mstock@umich.edu
 */

#include "vicmoc.h"
#ifdef _OPENMP
#include <omp.h>
#endif

// setup routines
int add_sharp_circular_blob (int,int,int,int,float**,float,float,float,float);
int add_smooth_circular_blob (int,int,int,int,float**,float,float,float,float);
int add_smooth_spherical_blob (int,int,int,int,int,int,float***,float,float,float,float,float);

// simulation routines
float step_forward_2d (int,int,int,int,int,int,int,int,float*,int,int,float***,float***,float***,int,float**,float,float*,int,float*,float,int,float,float,float***);

int create_baroclinic_vorticity_2d (int,int,int,int,float**,float**,float**,float,float*,float);
int create_boundary_vorticity_2d (int,int,int,int,float**,float**,float**,float);
float diffuse_scalar_2d (int,int,int,int,int,float**,float**,float,int,float**,float);
float variable_diffuse_scalar_2d (int,int,int,int,float**,float**,float**,int,float**,float);
int find_vels_2d (int,int,int,int,int,int,int,float*,float**,float**,float**,const int,float**,const float);
int find_gradient_of_scalar_2d (int,int,int,int,float**,float**,float,float**,float);
int find_gradient_of_scalar_2nd_2d (int,int,int,int,float**,float**,float**,float,float**,float);
int find_curl_of_vel_2d (int,int,int,int,float**,float**,float**);
int moc_advect_2d (int,int,int,int,float**,float**,float**,float*,float***,float***,float,int,int);
int find_open_boundary_psi (int,int,float**,float,float*,float**);
int find_biot_savart (int,float,float,int,int,float,float,float**,float*);
void recursive_biot_savart (const float, const float, const int, const int, const int, const float, const float, float***, float*);

// interpolation routines
float interpolate_array_using_CIC_2d (int,int,int,int,float***,float,float,int,float*);
int interpolate_vel_using_CIC_2d (int,int,int,int,float**,float**,float,float,float*,float*);

float interpolate_array_using_TSC_2d (int,int,int,int,float***,float,float,int,float*);
int interpolate_vel_using_TSC_2d (int,int,int,int,float**,float**,float,float,float*,float*);

float interpolate_array_using_M4p_2d (int,int,int,int,float**,float***,float,float,int,float*);
int interpolate_vel_using_M4p_2d (int,int,int,int,float**,float**,float**,float,float,float*,float*);


/*
 * For whatever reason, add a chunk of scalar stuff
 */
int add_sharp_circular_blob (int nx,int ny,int xbdry,int ybdry,float **addto,float xpos,float ypos,float rad,float val) {

   int ix,iy;
   float px,py;

   for (ix=0; ix<nx; ix++) {
      px = (float)ix/(float)(nx-1);
      for (iy=0; iy<ny; iy++) {
         py = (float)iy/(float)(nx-1);
         if (sqrt(pow(SCALE*(px-xpos),2)+pow(py-ypos,2)) < rad) addto[ix][iy] += val;
      }
   }

   return(0);
}


/*
 * For whatever reason, add a chunk of scalar stuff
 */
int add_smooth_circular_blob(int nx,int ny,int xbdry,int ybdry,float **addto,float xpos,float ypos,float rad,float val) {

   int ix,iy;
   float px,py;
   float rad_inner = 0.8*rad;
   float dist;

   for (ix=0; ix<nx; ix++) {
      px = (float)ix/(float)(nx-1);
      for (iy=0; iy<ny; iy++) {
         py = (float)iy/(float)(nx-1);
         dist = sqrt(pow(SCALE*(px-xpos),2)+pow(py-ypos,2));
         if (dist < rad_inner) {
            addto[ix][iy] += val;
         } else if (dist < rad) {
            addto[ix][iy] += val*(0.5+0.5*cos(M_PI*(dist-rad_inner)/(rad-rad_inner)));
         }
      }
   }

   return(0);
}


/*
 * For whatever reason, add a chunk of scalar stuff
 */
int add_smooth_spherical_blob(int nx,int ny,int nz,int xbdry,int ybdry,int zbdry,
   float ***addto,float xpos,float ypos,float zpos,float rad,float val) {

   int ix,iy,iz;
   float px,py,pz;
   float rad_inner = 0.8*rad;
   float dist;

   for (ix=0; ix<nx; ix++) {
      px = (float)ix/(float)(nx-1);
      for (iy=0; iy<ny; iy++) {
         py = (float)iy/(float)(nx-1);
         for (iz=0; iz<nz; iz++) {
            pz = (float)iz/(float)(nx-1);
            dist = sqrt(pow(px-xpos,2)+pow(py-ypos,2)+pow(pz-zpos,2));
            if (dist < rad_inner) {
               addto[ix][iy][iz] += val;
            } else if (dist < rad) {
               addto[ix][iy][iz] += val*(0.5+0.5*cos(M_PI*(dist-rad_inner)/(rad-rad_inner)));
            }
         }
      }
   }

   return(0);
}


/*
 * Procedure to take one step forward, given only vorticity
 * and scalar fraction
 */
float step_forward_2d (int silent, int step, int isStam, int mocOrder,
    int nx, int ny, int xbdry, int ybdry,
    float *freestream, int recalc_vel, int move_colors,
    float ***u, float ***a, float ***t,
    const int use_MASK, float **mask, const float maskerr,
    float sc_diffus[MAX_SCALARS],
    int gravtype, float *grav,
    float dt,
    int use_strong_strat, float bn, float dens_ratio, float ***acc) {

  int i,j,k;
  float Re,coeff;
  char outfileroot[70];

  // ----------------------------
  if (isStam) {

    // Use the method initially published in Jos Stam's "Stable Fluids"
    //  paper of 1998 whereupon the new velocity field is calculated by
    //  tracking cell centers one time step backwards and reading the
    //  interpolated velocity at that location.

    // There is no need to create boundary vorticity, as the no-slip
    //  condition is enforced during the velocity solution

    // But, how do we enforce baroclinicity?

    // regardless of the above, copy 'u' to 'a' so that we can diffuse it later
    (void) copy_2d_field(u[XV], a[XV], nx, ny);
    (void) copy_2d_field(u[YV], a[YV], nx, ny);

  } else {

    // This uses the new, unpublished method by Mark Stock whereupon the
    //  new *vorticity* field is updated in the same manner as Stam's
    //  velocity field.

    // create new vorticity based on density differences (gradient of scalar)
    if (a[SF] != NULL)
      create_baroclinic_vorticity_2d (nx,ny,xbdry,ybdry,a[SF],mask,a[W2],dt,grav,bn);

    // create new vorticity due to viscous boundaries, etc,
    // sets 'a' vorticity at boundaries
    create_boundary_vorticity_2d (nx,ny,xbdry,ybdry,u[XV],u[YV],a[W2],dt);

  }

  // diffuse each scalar individually (a[] becomes t[])
  for (k=0; k<MAX_SCALARS; k++) {
    // first, do the arrays exist?
    if (a[k] != NULL && t[k] != NULL) {
      // next, do we even diffuse them?
      if (TRUE) {
        // finally, is the diffusivity uniform or variable?
        if (sc_diffus[k] < 0.0) {
          // run the variable viscosity routine (always use MD for now)
          coeff = variable_diffuse_scalar_2d (nx,ny,xbdry,ybdry,a[k],t[k],
                                     a[MD],use_MASK,mask,dt);
        } else {
          // run the uniform viscosity routine
          coeff = diffuse_scalar_2d (nx,ny,xbdry,ybdry,k,a[k],t[k],
                                     sc_diffus[k],use_MASK,mask,dt);
        }
      } else {
        // if we don't diffuse it, we still need to copy the values
        for (i=0;i<nx;i++) for (j=0;j<ny;j++) {
          t[k][i][j] = a[k][i][j];
        }
        coeff = 1.0;
      }
    }
  }

  // again, split on method
  if (isStam) {

    // must copy 't' back to 'u' so that it can be used in moc advection
    (void) copy_2d_field(t[XV], u[XV], nx, ny);
    (void) copy_2d_field(t[YV], u[YV], nx, ny);

    //sprintf(outfileroot,"u_%04d",step+30);
    //write_output(outfileroot,nx,ny,u[XV],-0.25,0.5,1);
    //sprintf(outfileroot,"v_%04d",step+30);
    //write_output(outfileroot,nx,ny,u[YV],-0.25,0.5,1);

  }

  // project forward to find the new fields (both)
  if (!silent) fprintf(stderr,"  now in moc_advect_2d\n"); fflush(stderr);
  moc_advect_2d (nx,ny,xbdry,ybdry,mask,u[XV],u[YV],freestream,t,a,dt,mocOrder,move_colors);

  // split on method for advection step
  if (isStam) {

    // apply no-slip wall boundary conditions
    if (ybdry == WALL) {
      j = ny-1;
      for (i=0;i<nx;i++) {
         a[XV][i][j] = 0.;
         a[YV][i][j] = 0.;
      }
      j = 0;
      for (i=0;i<nx;i++) {
         a[XV][i][j] = 0.;
         a[YV][i][j] = 0.;
      }
    }
    if (xbdry == WALL) {
      i = 0;
      for (j=0;j<ny;j++) {
         a[XV][i][j] = 0.;
         a[YV][i][j] = 0.;
      }
      i = nx-1;
      for (j=0;j<ny;j++) {
         a[XV][i][j] = 0.;
         a[YV][i][j] = 0.;
      }
    }

    // already used 'u' to find 'a' field from 't' field

    // calculate the vorticity field from 'a' field
    find_curl_of_vel_2d (nx,ny,xbdry,ybdry,a[XV],a[YV],a[W2]);

    // project forward to find the new velocity field
    // project 'a' velocity field onto a divergence-free vel field 'u'

    // make sure that the true velocities are in 'u' before going back

  } else {

    // create new vorticity due to viscous boundaries, etc,
    // sets 'a' vorticity at boundaries
    create_boundary_vorticity_2d (nx,ny,xbdry,ybdry,u[XV],u[YV],a[W2],dt);

  }

  // find the velocity field from the vorticity field
  if (recalc_vel) {
    if (!silent) fprintf(stderr,"  now in find_vels_2d\n"); fflush(stderr);
    find_vels_2d (silent,step,isStam,nx,ny,xbdry,ybdry,freestream,u[XV],u[YV],a[W2],use_MASK,mask,maskerr);
  }

  // calculate effective Reynolds number
  Re = pow((float)(nx-1),2)*dt/coeff;

  // return Re, why? I dunno.
  return(Re);
}


/*
 * create vorticity at areas of scalar gradient
 */
int create_baroclinic_vorticity_2d (int nx,int ny,int xbdry,int ybdry,
      float **scalar,float **mask,float **vort,float dt,float *g,float bn) {

   int i,j;
   //float g[2] = {0.0,-1.0};
   static float **grad[2];
   static int set_grad = FALSE;

   // fprintf(stderr,"in create_baroclinic_vorticity_2d\n"); fflush(stderr);

   // if no boussinesq number (bn ~ 0), skip this!
   if (fabs(bn) < EPSILON) return(0);

   // allocate memory for the gradient
   if (!set_grad) {
      grad[0] = allocate_2d_array_f(nx,ny);
      grad[1] = allocate_2d_array_f(nx,ny);
      set_grad = TRUE;
   }

   // first, find the scalar gradient
   find_gradient_of_scalar_2nd_2d (nx,ny,xbdry,ybdry,scalar,mask,grad[0],bn,grad[1],bn);

   // then, cross with the gravity vector
   #pragma omp parallel for private(i,j)
   for (i=0;i<nx;i++) {
      for (j=0;j<ny;j++) {
         vort[i][j] += dt*(grad[0][i][j]*g[1] - grad[1][i][j]*g[0]);
      }
   }

   // nah, just keep it. it's OK. it needs to diffuse, anyway
   return(0);
}


/* 
 * create vorticity at solid boundaries - 2D
 */
int create_boundary_vorticity_2d (int nx,int ny,int xbdry,int ybdry,
      float **u,float **v,float **vort,float dt) {

   int skip_it = FALSE;
   int i,j;
   int nxm1 = nx-1;
   int nym1 = ny-1;
   int nxm2 = nx-2;
   int nym2 = ny-2;
   int nxm3 = nx-3;
   int nym3 = ny-3;
   float hxi,hyi;
   //float mult = 1.0*nxm1*dt;	// 1.0 seems to be correct

   // old way (not right)
/*
   hxi = (float)(nxm1);
   hyi = (float)(nxm1);
   if (xbdry == WALL) {
      for (j=0;j<ny;j++) {
         vort[0][j] += mult*v[0][j];
         vort[nxm1][j] -= mult*v[nxm1][j];
      }
   }
   if (ybdry == WALL) {
      for (i=0;i<nx;i++) {
         vort[i][0] -= mult*u[i][0];
         vort[i][nym1] += mult*u[i][nym1];
      }
   }
*/

   if (skip_it) {

   if (xbdry == WALL) {
      for (j=0; j<ny; j++) {
         vort[0][j] = 0.;
         vort[nxm1][j] = 0.;
      }
   }
   if (ybdry == WALL) {
      for (i=0; i<nx; i++) {
         vort[i][0] = 0.;
         vort[i][nym1] = 0.;
      }
   }

   } else {

   // new way, from CottetVM - OOOOooooh, I was using slip walls!
   // hxi = dt*0.5*(float)(nxm1);
   // hyi = dt*0.5*(float)(nxm1);
   //hxi = dt*0.5*(float)(nxm1);
   //hyi = -dt*0.5*(float)(nxm1);
   // DUH! When I redefine vorticity (as I need to do to force the wall
   //  velocity to zero) I don't need to multiply by dt! Dummy.
   //hxi = 0.5*(float)(nxm1)/SCALE;
   //hyi = -0.5*(float)(nxm1);
   if (nx > ny) {
      hxi = 0.5*(float)(nxm1);
      hyi = -0.5*(float)(nxm1);
   } else {
      hxi = 0.5*(float)(nym1);
      hyi = -0.5*(float)(nym1);
   }
   //hxi = dt*0.5*(float)(nxm1)/SCALE;
   //hyi = dt*-0.5*(float)(nxm1);
   if (xbdry == WALL) {
      for (j=0; j<ny; j++) {
         // vort[0][j] = (-1.*v[2][j]+4.*v[1][j]-3.*v[0][j])*hxi;
         // vort[nxm1][j] = (v[nxm3][j]-4.*v[nxm2][j]+3.*v[nxm1][j])*hxi;
         vort[0][j] = (-1.*v[2][j]+4.*v[1][j])*hxi;
         vort[nxm1][j] = (v[nxm3][j]-4.*v[nxm2][j])*hxi;
         // vort[0][j] += (-1.*v[2][j]+4.*v[1][j])*hxi;
         // vort[nxm1][j] += (v[nxm3][j]-4.*v[nxm2][j])*hxi;
         // vort[0][j] += (-1.*v[2][j]+4.*v[1][j]-3.*v[0][j])*hxi;
         // vort[nxm1][j] += (v[nxm3][j]-4.*v[nxm2][j]+3.*v[nxm1][j])*hxi;
      }
   }
   if (ybdry == WALL) {
      for (i=0; i<nx; i++) {
         // vort[i][0] = (-1.*u[i][2]+4.*u[i][1]-3.*u[i][0])*hyi;
         // vort[i][nym1] = (u[i][nym3]-4.*u[i][nym2]+3.*u[i][nym1])*hyi;
         vort[i][0] = (-1.*u[i][2]+4.*u[i][1])*hyi;
         vort[i][nym1] = (u[i][nym3]-4.*u[i][nym2])*hyi;
         // vort[i][0] += (-1.*u[i][2]+4.*u[i][1])*hyi;
         // vort[i][nym1] += (u[i][nym3]-4.*u[i][nym2])*hyi;
         // vort[i][0] += (-1.*u[i][2]+4.*u[i][1]-3.*u[i][0])*hyi;
         // vort[i][nym1] += (u[i][nym3]-4.*u[i][nym2]+3.*u[i][nym1])*hyi;
      }
   }

   }

   return(0);
}


/* 
 * diffuse scalar field according to laplacian of vorticity, return
 * to a different scalar field
 *
 * diffus is the scalar diffusivity, like nu for vorticity: 1/Re
 */
float diffuse_scalar_2d (int nx,int ny,int xbdry,int ybdry,
      int vartype,float **in,float **out,
      float diffus,int use_MASK,float **mask,float dt) {

   int i,j;
   int nxm1 = nx-1;
   int nym1 = ny-1;
   int istep,numsteps;
   float multx,multy;
   float mult,omm;

   if (nx > ny) {
      multx = dt*diffus*pow(nxm1,2);
      multy = dt*diffus*pow(nxm1,2);
   } else {
      multx = dt*diffus*pow(nym1,2);
      multy = dt*diffus*pow(nym1,2);
   }

   // if mult>0.125 (stability requirement) then run multiple substeps
   //   each with smaller dt
   if (multx > multy) mult = multx;
   else mult = multy;
   numsteps = 1 + (int)(mult/0.125);
   multx = multx/(float)(numsteps);
   multy = multy/(float)(numsteps);

   for (istep=0; istep<numsteps; istep++) {

   // do the middle part of the field, regardless of the BCs
   if (use_MASK) {
     if (vartype == W2) {
       // this is vorticity, allow it to diffuse from masks
       #pragma omp parallel for private(i,j)
       for (i=1;i<nxm1;i++) {
         for (j=1;j<nym1;j++) {
           // diffuse vort from masks, but not into masks
           out[i][j] = in[i][j]
                     + mask[i][j]*multx*(in[i+1][j]+in[i-1][j]-2.*in[i][j])
                     + mask[i][j]*multy*(in[i][j+1]+in[i][j-1]-2.*in[i][j]);
         }
       }
     } else {
       // this is not vorticity, prevent it from diffusing from masks
       #pragma omp parallel for private(i,j)
       for (i=1;i<nxm1;i++) {
         for (j=1;j<nym1;j++) {
           // avoid diffusing anything from masks
           // first attempt 2014-10-20
           //out[i][j] = in[i][j]
           //          + mask[i+1][j]*mask[i][j]*mask[i-1][j]*multx*(in[i+1][j]+in[i-1][j]-2.*in[i][j])
           //          + mask[i][j+1]*mask[i][j]*mask[i][j-1]*multy*(in[i][j+1]+in[i][j-1]-2.*in[i][j]);
           // see notes 2014-10-21
           out[i][j] = in[i][j]
                     + mask[i][j]*multx*(mask[i+1][j]*(in[i+1][j]-in[i][j]) + mask[i-1][j]*(in[i-1][j]-in[i][j]))
                     + mask[i][j]*multy*(mask[i][j+1]*(in[i][j+1]-in[i][j]) + mask[i][j-1]*(in[i][j-1]-in[i][j]));
         }
       }
     }
   } else {
     #pragma omp parallel for private(i,j)
     for (i=1;i<nxm1;i++) {
       for (j=1;j<nym1;j++) {
         out[i][j] = in[i][j]
                   + multx*(in[i+1][j]+in[i-1][j]-2.*in[i][j])
                   + multy*(in[i][j+1]+in[i][j-1]-2.*in[i][j]);
       }
     }
   }

   // walls, 2nd order Laplacian
   // then solve for the wall and periodic bdry
   // new: do not diffuse on OPEN boundaries, because flow may be coming in
   if (xbdry == WALL) {
      for (j=1;j<nym1;j++) {
         out[0][j] = in[0][j]
                   + multx*(-in[2][j]+6.*in[1][j]-5.*in[0][j])
                   + multy*(in[0][j+1]+in[0][j-1]-2.*in[0][j]);
         out[nxm1][j] = in[nxm1][j]
                   + multx*(-in[nx-3][j]+6.*in[nx-2][j]-5.*in[nxm1][j])
                   + multy*(in[nxm1][j+1]+in[nxm1][j-1]-2.*in[nxm1][j]);
      }
   } else if (xbdry == PERIODIC) {
      for (j=1;j<nym1;j++) {
         out[0][j] = in[0][j]
                   + multx*(in[1][j]+in[nx-2][j]-2.*in[0][j])
                   + multy*(in[0][j+1]+in[0][j-1]-2.*in[0][j]);
         out[nxm1][j] = out[0][j];
      }
   } else {	// OPEN
      for (j=1;j<nym1;j++) {
         out[0][j] = in[0][j];
         out[nxm1][j] = in[nxm1][j];
      }
   }
   if (ybdry == WALL) {
      for (i=1;i<nxm1;i++) {
         out[i][0] = in[i][0]
                   + multx*(in[i+1][0]+in[i-1][0]-2.*in[i][0])
                   + multy*(-in[i][2]+6.*in[i][1]-5.*in[i][0]);
         out[i][nym1] = in[i][nym1]
                   + multx*(in[i+1][nym1]+in[i-1][nym1]-2.*in[i][nym1])
                   + multy*(-in[i][ny-3]+6.*in[i][ny-2]-5.*in[i][nym1]);
      }
   } else if (ybdry == PERIODIC) {
      for (i=1;i<nxm1;i++) {
         out[i][0] = in[i][0]
                   + multx*(in[i+1][0]+in[i-1][0]-2.*in[i][0])
                   + multy*(in[i][1]+in[i][ny-2]-2.*in[i][0]);
         out[i][nym1] = out[i][0];
      }
   } else {	// OPEN
      for (i=1;i<nxm1;i++) {
         out[i][0] = in[i][0];
         out[i][nym1] = in[i][nym1];
      }
   }

   // corners, 2nd order Laplacian
   if (xbdry == WALL && ybdry == WALL) {
      out[0][0] = in[0][0]
                   + multx*(-in[2][0]+6.*in[1][0]-5.*in[0][0])
                   + multy*(-in[0][2]+6.*in[0][1]-5.*in[0][0]);
      out[0][nym1] = in[0][nym1]
                   + multx*(-in[2][nym1]+6.*in[1][nym1]-5.*in[0][nym1])
                   + multy*(-in[0][ny-3]+6.*in[0][ny-2]-5.*in[0][nym1]);
      out[nxm1][0] = in[nxm1][0]
                   + multx*(-in[nx-3][0]+6.*in[nx-2][0]-5.*in[nxm1][0])
                   + multy*(-in[nxm1][2]+6.*in[nxm1][1]-5.*in[nxm1][0]);
      out[nxm1][nym1] = in[nxm1][nym1]
                   + multx*(-in[nx-3][nym1]+6.*in[nx-2][nym1]-5.*in[nxm1][nym1])
                   + multy*(-in[nxm1][ny-3]+6.*in[nxm1][ny-2]-5.*in[nxm1][nym1]);
   } else if (xbdry == WALL && ybdry == PERIODIC) {
      // t-kernel 2nd derivatives for two, then copy
      out[0][0] = in[0][0]
                   + multx*(-in[2][0]+6.*in[1][0]-5.*in[0][0])
                   + multy*(in[0][1]+in[0][ny-2]-2.*in[0][0]);
      out[nxm1][0] = in[nxm1][0]
                   + multx*(-in[nx-3][0]+6.*in[nx-2][0]-5.*in[nxm1][0])
                   + multy*(in[nxm1][1]+in[nxm1][ny-2]-2.*in[nxm1][0]);
      out[0][nym1] = out[0][0];
      out[nxm1][nym1] = out[nxm1][0];
   } else if (xbdry == PERIODIC && ybdry == WALL) {
      // t-kernel 2nd derivatives for two, then copy
      out[0][0] = in[0][0]
                   + multx*(in[1][0]+in[nx-2][0]-2.*in[0][0])
                   + multy*(-in[0][2]+6.*in[0][1]-5.*in[0][0]);
      out[0][nym1] = in[0][nym1]
                   + multx*(in[1][nym1]+in[nx-2][nym1]-2.*in[0][nym1])
                   + multy*(-in[0][ny-3]+6.*in[0][ny-2]-5.*in[0][nym1]);
      out[nxm1][0] = out[0][0];
      out[nxm1][nym1] = out[0][nym1];
   } else if (xbdry == PERIODIC && ybdry == PERIODIC) {
      out[0][0] = in[0][0]
                   + multx*(in[1][0]+in[nx-2][0]-2.*in[0][0])
                   + multy*(in[0][1]+in[0][ny-2]-2.*in[0][0]);
      out[0][nym1] = out[0][0];
      out[nxm1][0] = out[0][0];
      out[nxm1][nym1] = out[0][0];
   } else {	// all OPEN
      out[0][0] = in[0][0];
      out[0][nym1] = in[0][nym1];
      out[nxm1][0] = in[nxm1][0];
      out[nxm1][nym1] = in[nxm1][nym1];
   }

   // and, for subsequent steps to act on the right variables, swap
   for (i=0;i<nx;i++) {
      for (j=0;j<ny;j++) {
         in[i][j] = out[i][j];
      }
   }

   }	// end loop over numsteps

   return(mult*numsteps);
}


/* 
 * diffuse scalar field according to laplacian of vorticity, return
 * to a different scalar field
 *
 * where mask is 1.0, do not diffuse, where mask is 0.0, fully diffuse
 *
 * diffus is the scalar diffusivity, like nu for vorticity: 1/Re
 */
float variable_diffuse_scalar_2d (int nx,int ny,int xbdry,int ybdry,
      float **in,float **out,
      float **mu,int use_MASK,float **mask,float dt) {

   int i,j;
   int nxm1 = nx-1;
   int nym1 = ny-1;
   int istep,numsteps;
   float diffus;
   float mult,this,omm;

   // find maximum diffusivity
   diffus = 0.0;
   for (i=0;i<nx;i++) {
      for (j=0;j<ny;j++) {
         if (mu[i][j] > diffus) diffus = mu[i][j];
      }
   }

   // test multiple in both directions
   if (nx > ny) {
      mult = dt*diffus*pow(nxm1,2);
   } else {
      mult = dt*diffus*pow(nym1,2);
   }

   // if mult>0.125 (stability requirement) then run multiple substeps
   //   each with smaller dt
   numsteps = 1 + (int)(mult/0.125);
   // also divide by diffus (max mu) because we multiply by mu later
   mult = mult/(diffus*numsteps);
   if (numsteps > 1) {
      fprintf(stdout,"  diffusing %d steps\n",numsteps);
      fflush(stdout);
   }

   for (istep=0; istep<numsteps; istep++) {

   // do the middle part of the field, regardless of the BCs
   //if (use_MASK) {
   if (FALSE) {
     #pragma omp parallel for private(i,j)
     for (i=1;i<nxm1;i++) {
       for (j=1;j<nym1;j++) {
         omm = 1.-mask[i][j];
         this = mult*mu[i][j];
         out[i][j] = in[i][j] + omm*this*(in[i+1][j]+in[i-1][j]-2.*in[i][j]) + omm*this*(in[i][j+1]+in[i][j-1]-2.*in[i][j]);
       }
     }
   } else {
     #pragma omp parallel for private(i,j)
     for (i=1;i<nxm1;i++) {
       for (j=1;j<nym1;j++) {
         this = mult*mu[i][j];
         out[i][j] = in[i][j] + this*(in[i+1][j]+in[i-1][j]-2.*in[i][j]) + this*(in[i][j+1]+in[i][j-1]-2.*in[i][j]);
       }
     }
   }

   // then solve for the wall and periodic bdry
   if (xbdry == WALL || xbdry == OPEN) {
      // for (j=0;j<ny;j++) {
      for (j=1;j<nym1;j++) {
         // out[0][j] = 0.0;
         //out[0][j] = in[0][j] + mult*(in[0][j+1]+in[0][j-1]-in[2][j]+6*in[1][j]-7*in[0][j]);
         this = mult*mu[0][j];
         out[0][j] = in[0][j] + this*(-in[2][j]+6.*in[1][j]-5.*in[0][j]) + this*(in[0][j+1]+in[0][j-1]-2.*in[0][j]);
         // out[nxm1][j] = 0.0;
         //out[nxm1][j] = in[nxm1][j] + mult*(in[nxm1][j+1]+in[nxm1][j-1]-in[nx-3][j]+6*in[nx-2][j]-7*in[nxm1][j]);
         this = mult*mu[nxm1][j];
         out[nxm1][j] = in[nxm1][j] + this*(-in[nx-3][j]+6.*in[nx-2][j]-5.*in[nxm1][j]) + this*(in[nxm1][j+1]+in[nxm1][j-1]-2.*in[nxm1][j]);
      }
   } else if (xbdry == PERIODIC) {
      for (j=1;j<nym1;j++) {
         //out[0][j] = in[0][j] + mult*(in[1][j]+in[nx-2][j]+in[0][j+1]+in[0][j-1]-4.0*in[0][j]);
         this = mult*mu[0][j];
         out[0][j] = in[0][j] + this*(in[1][j]+in[nx-2][j]-2.*in[0][j]) + this*(in[0][j+1]+in[0][j-1]-2.*in[0][j]);
         out[nxm1][j] = out[0][j];
      }
   }
   if (ybdry == WALL || ybdry == OPEN) {
      // for (i=0;i<nx;i++) {
      for (i=1;i<nxm1;i++) {
         // out[i][0] = 0.0;
         //out[i][0] = in[i][0] + mult*(in[i+1][0]+in[i-1][0]-in[i][2]+6*in[i][1]-7*in[i][0]);
         this = mult*mu[i][0];
         out[i][0] = in[i][0] + this*(in[i+1][0]+in[i-1][0]-2.*in[i][0]) + this*(-in[i][2]+6.*in[i][1]-5.*in[i][0]);
         // out[i][nym1] = 0.0;
         //out[i][nym1] = in[i][nym1] + mult*(in[i+1][nym1]+in[i-1][nym1]-in[i][ny-3]+6*in[i][ny-2]-7*in[i][nym1]);
         this = mult*mu[i][nym1];
         out[i][nym1] = in[i][nym1] + this*(in[i+1][nym1]+in[i-1][nym1]-2.*in[i][nym1]) + this*(-in[i][ny-3]+6.*in[i][ny-2]-5.*in[i][nym1]);
      }
   } else if (ybdry == PERIODIC) {
      for (i=1;i<nxm1;i++) {
         //out[i][0] = in[i][0] + mult*(in[i+1][0]+in[i-1][0]+in[i][1]+in[i][ny-2]-4.0*in[i][0]);
         this = mult*mu[i][0];
         out[i][0] = in[i][0] + this*(in[i+1][0]+in[i-1][0]-2.*in[i][0]) + this*(in[i][1]+in[i][ny-2]-2.*in[i][0]);
         out[i][nym1] = out[i][0];
      }
   }
   if ((xbdry == WALL && ybdry == WALL) || (xbdry == OPEN && ybdry == OPEN)) {
      // corner 2nd derivatives for all 4
      //out[0][0] = in[0][0] + mult*(-in[2][0]-in[0][2]+6*in[1][0]+6*in[0][1]-10*in[0][0]);
      this = mult*mu[0][0];
      out[0][0] = in[0][0] + this*(-in[2][0]+6.*in[1][0]-5.*in[0][0]) + this*(-in[0][2]+6.*in[0][1]-5.*in[0][0]);
      //out[0][nym1] = in[0][nym1] + mult*(-in[2][nym1]-in[0][ny-3]+6*in[1][nym1]+6*in[0][ny-2]-10*in[0][nym1]);
      this = mult*mu[0][nym1];
      out[0][nym1] = in[0][nym1] + this*(-in[2][nym1]+6.*in[1][nym1]-5.*in[0][nym1]) + this*(-in[0][ny-3]+6.*in[0][ny-2]-5.*in[0][nym1]);
      //out[nxm1][0] = in[nxm1][0] + mult*(-in[nx-3][0]-in[nxm1][2]+6*in[nx-2][0]+6*in[nxm1][1]-10*in[nxm1][0]);
      this = mult*mu[nxm1][0];
      out[nxm1][0] = in[nxm1][0] + this*(-in[nx-3][0]+6.*in[nx-2][0]-5.*in[nxm1][0]) + this*(-in[nxm1][2]+6.*in[nxm1][1]-5.*in[nxm1][0]);
      //out[nxm1][nym1] = in[nxm1][nym1] + mult*(-in[nx-3][nym1]-in[nxm1][ny-3]+6*in[nx-2][nym1]+6*in[nxm1][ny-2]-10*in[nxm1][nym1]);
      this = mult*mu[nxm1][nym1];
      out[nxm1][nym1] = in[nxm1][nym1] + this*(-in[nx-3][nym1]+6.*in[nx-2][nym1]-5.*in[nxm1][nym1]) + this*(-in[nxm1][ny-3]+6.*in[nxm1][ny-2]-5.*in[nxm1][nym1]);
   } else if (xbdry == WALL && ybdry == PERIODIC) {
      // t-kernel 2nd derivatives for two, then copy
      //out[0][0] = in[0][0] + mult*(in[0][1]+in[0][ny-2]-in[2][0]+6*in[1][0]-7*in[0][0]);
      this = mult*mu[0][0];
      out[0][0] = in[0][0] + this*(-in[2][0]+6.*in[1][0]-5.*in[0][0]) + this*(in[0][1]+in[0][ny-2]-2.*in[0][0]);
      //out[nxm1][0] = in[nxm1][0] + mult*(in[nxm1][1]+in[nxm1][ny-2]-in[nx-3][0]+6*in[nx-2][0]-7*in[nxm1][0]);
      this = mult*mu[nxm1][0];
      out[nxm1][0] = in[nxm1][0] + this*(-in[nx-3][0]+6.*in[nx-2][0]-5.*in[nxm1][0]) + this*(in[nxm1][1]+in[nxm1][ny-2]-2.*in[nxm1][0]);
      out[0][nym1] = out[0][0];
      out[nxm1][nym1] = out[nxm1][0];
   } else if (xbdry == PERIODIC && ybdry == WALL) {
      // t-kernel 2nd derivatives for two, then copy
      //out[0][0] = in[0][0] + mult*(in[1][0]+in[nx-2][0]-in[0][2]+6*in[0][1]-7*in[0][0]);
      this = mult*mu[0][0];
      out[0][0] = in[0][0] + this*(in[1][0]+in[nx-2][0]-2.*in[0][0]) + this*(-in[0][2]+6.*in[0][1]-5.*in[0][0]);
      //out[0][nym1] = in[0][nym1] + mult*(in[1][nym1]+in[nx-2][nym1]-in[0][ny-3]+6*in[0][ny-2]-7*in[0][nym1]);
      this = mult*mu[0][nym1];
      out[0][nym1] = in[0][nym1] + this*(in[1][nym1]+in[nx-2][nym1]-2.*in[0][nym1]) + this*(-in[0][ny-3]+6.*in[0][ny-2]-5.*in[0][nym1]);
      out[nxm1][0] = out[0][0];
      out[nxm1][nym1] = out[0][nym1];
   } else if (xbdry == PERIODIC && ybdry == PERIODIC) {
      //out[0][0] = in[0][0] + mult*(in[0][1]+in[0][ny-2]+in[1][0]+in[nx-2][0]-4*in[0][0]);
      this = mult*mu[0][0];
      out[0][0] = in[0][0] + this*(in[1][0]+in[nx-2][0]-2.*in[0][0]) + this*(in[0][1]+in[0][ny-2]-2.*in[0][0]);
      out[0][nym1] = out[0][0];
      out[nxm1][0] = out[0][0];
      out[nxm1][nym1] = out[0][0];
   } else {
      fprintf(stderr,"ERROR (variable_diffuse_scalar_2d): should not get here\n");
      exit(0);
   }

   // and, for subsequent steps to act on the right variables, swap
   for (i=0;i<nx;i++) {
      for (j=0;j<ny;j++) {
         in[i][j] = out[i][j];
      }
   }

   //if (istep > 1) fprintf(stdout,"."); fflush(stdout);
   }	// end loop over numsteps
   //if (numsteps > 1) fprintf(stdout,"\n"); fflush(stdout);

   return(mult*diffus*numsteps);
}


/*
 * find_vels_2d takes a vorticity field, and solves for the velocity field
 *
 * mask is a flow mask, where 1.0 means solid and 0.0 means open
 *
 * if driven cavity problem, force the lid to vel=1
 */
int find_vels_2d (int silent, int step,const int isStam,const int nx,const int ny,
      const int xbdry,const int ybdry,float *freestream,
      float **u,float **v,float **vort,
      const int use_MASK,float **mask,const float maskerr) {

   int use_multigrid = TRUE;	// use MUDPACK solver mud2sp.f
   //int use_multigrid = FALSE;	// use FISHPAK solver gr2.c
   int driven_cavity = FALSE;
   int i,j;
   // these are for hwscrt, FISHPAK
   int nxm1 = nx-1;
   int nym1 = ny-1;
   int nxp = nx;
   int nyp = ny;
   int bdcx,bdcy;
   float bdvxs[ny];
   float bdvxf[ny];
   float bdvys[nx];
   float bdvyf[nx];
   float elmbda = 0.0;
   float pert;
   // these are for mud2sp, MUDPACK
   static int iparm[17];
   static float fparm[6];
   static char bndyc[10];
   static char cofx[10];
   static char cofy[10];
   static float **rhs;
   static int mgopt[4];
   static int must_initialize = TRUE;
   float muderr = 1.e-4;
   // these are for both
   int ierr;
   long int iworksize = 3.67*nx*ny+100*nx+2000;	// OK for all?
   static float *work;	// increasing this does not fix hwscrt's problem
   float xs,xf,ys,yf;
   static float **psi;
   static int allocated_arrays = FALSE;
   int istep,maxstep;
   float maxvort,temp,maxcorr;
   static float **cu;
   static float **cv;
   static float **cw;
   static float **cw_old;
   char outfileroot[150];

   //fprintf(stderr,"in find_vels_2d\n"); fflush(stderr);

   // allocate memory for the streamfunction, and the boundary
   //    vorticity arrays
   if (!allocated_arrays) {
      psi = allocate_2d_array_f(nx,ny);
      work = allocate_1d_array_f(iworksize);

      // and allocate rhs, if using multigrid solver
      if (use_multigrid)
         rhs = allocate_2d_array_f(nx,ny);

      // allocate arrays for correction velocities
      if (use_MASK) {
         cu = allocate_2d_array_f(nx,ny);
         cv = allocate_2d_array_f(nx,ny);
         cw = allocate_2d_array_f(nx,ny);
         cw_old = allocate_2d_array_f(nx,ny);
      }

      allocated_arrays = TRUE;
   }

   // set domain dimensions
   // how does SCALE work into this?
   xs = 0.0;
   ys = 0.0;
   if (nx > ny) {
     xf = 1.0;
     yf = ((float)(ny-1)/(float)(nx-1));
   } else {
     xf = ((float)(nx-1)/(float)(ny-1));
     yf = 1.0;
   }

   // set maximum number of iterations
   if (use_MASK) {
      maxstep = 1000;
   } else {
      maxstep = 1;
   }

   // begin iterative solution
   for (istep=0; istep<maxstep; istep++) {


   // split on solution technique

   // Use mudpack multigrid solver ----------------------------------------
   if (use_multigrid) {

      // only do this once
      if (must_initialize) {

         // Note that I swap x and y for everything relating to mud2sp
         // This is because 2D array ordering is opposite between C and F90

         // set the input parameters
         iparm[0] = 0;	// initialize the discretization

         // set the boundary conditions
         if (ybdry == PERIODIC) {
           iparm[1] = 0;
           iparm[2] = 0;
         } else if (ybdry == WALL || ybdry == OPEN) {
           iparm[1] = 1;
           iparm[2] = 1;
         }
         if (xbdry == PERIODIC) {
           iparm[3] = 0;
           iparm[4] = 0;
         } else if (xbdry == WALL || xbdry == OPEN) {
           iparm[3] = 1;
           iparm[4] = 1;
         }

         // find the proper integer dimensions
         // do x first
         //fprintf(stdout," nx is %d\n",nx);
         iparm[9] = ny;
         for (i=2; i<101; i++) {
           for (j=1; j<20; j++) {
             if (i*(pow(2,j))+1 == ny) {
               fprintf(stdout," ny (%d) is 1 + %d * 2^%d\n",ny,i,j);
               iparm[5] = i;
               iparm[7] = j+1;
               j = 100;
               i = 100;
             }
           }
         }
         if (i == 11) {
           fprintf(stderr,"Try an ny that is divisible by an integer from");
           fprintf(stderr," 2 to 10\n");
           fprintf(stderr," ny = %d\n",ny);
           fprintf(stderr,"Quitting.\n");
           exit(0);
         }
         if (j == 20) {
           fprintf(stderr,"ny = %d ?!? Are you sure?\n",nx);
           fprintf(stderr,"That's very large. 10k x 10k needs 10GB of RAM.\n");
           fprintf(stderr,"Quitting.\n");
           exit(0);
         }

         // do y next
         iparm[10] = nx;
         for (i=2; i<101; i++) {
           for (j=1; j<20; j++) {
             if (i*(pow(2,j))+1 == nx) {
               fprintf(stdout," nx (%d) is 1 + %d * 2^%d\n",nx,i,j);
               iparm[6] = i;
               iparm[8] = j+1;
               j = 100;
               i = 100;
             }
           }
         }
         if (i == 11) {
           fprintf(stderr,"Try an nx that is divisible by an integer from");
           fprintf(stderr," 2 to 10\n");
           fprintf(stderr," nx = %d\n",nx);
           fprintf(stderr,"Quitting.\n");
           exit(0);
         }
         if (j == 20) {
           fprintf(stderr,"nx = %d ?!? Are you sure?\n",ny);
           fprintf(stderr,"That's very large. 10k x 10k needs 10GB of RAM.\n");
           fprintf(stderr,"Quitting.\n");
           exit(0);
         }

         iparm[11] = 0;	// no initial guess for psi is provided
			// why can't we provide psi from the last step?
         iparm[12] = 100;	// max # multigrid cycles at finest res
         iparm[13] = 0;	// method of relaxation (0=point, 3=line x and y)
         iparm[14] = iworksize;	// size of workspace

         fparm[0] = ys;	// x start
         fparm[1] = yf;	// x end
         fparm[2] = xs;	// y start
         fparm[3] = xf;	// y end
         fparm[4] = muderr;	// tolerance criterion
         if (!silent) fprintf(stderr,"domain is %g < x < %g\n",xs,xf);
         if (!silent) fprintf(stderr,"          %g < y < %g\n",ys,yf);

         // always run a specific number of cycles
         //fparm[4] = 0.0;
         //iparm[12] = 10;

         mgopt[0] = 0;	// use default multigrid options

         // the boundary condition subroutine names
         // Apparently, the F77 solver calls these as "external"
         // I have no clue how to actually get them to work
         // For now, I have just replaced the calls in mud2sp_full.f with val=1.
         //   (look for "NOTE" or "cfx")
         (void) strcpy(cofx,"cofx");
         (void) strcpy(cofy,"cofy");
         // bndyc only used for mixed BCs - test Dirichlet BC first
         (void) strcpy(bndyc,"bndyc");

         // this flags the routine to discretize the problem and check input
         if (!silent) fprintf(stdout,"Running initialization mud2sp\n");
         mud2sp_(iparm, fparm, work, cofx, cofy,
                 bndyc, rhs[0], psi[0], mgopt, &ierr);

         // catch an error in the setup or input
         if (ierr != 0) {
            fprintf(stderr,"ERROR (mud2sp_): ierr = %d\n",ierr);
            if (ierr == 9) fprintf(stderr,"  iparm[14] (%d) too small, must be > %d\n",iparm[14],iparm[15]);
            fprintf(stderr,"Quitting.\n");
            exit(0);
         }

         // read workspace requirements!
         if (!silent) fprintf(stderr,"work array needs %d\n",iparm[15]);
         if (!silent) fprintf(stderr,"  we currently allocate %d\n",iworksize);
         if (iparm[15] > iworksize) {
           free_1d_array_f(work);
           work = allocate_1d_array_f(iworksize);
         }

         // don't do this again
         must_initialize = FALSE;

         // rest the intl parameter
         iparm[0] = 1;

      } else {	// not must_initialize

         // set the input parameters
         iparm[0] = 1;	// run a regular solution
      }


      // now call the subroutine for real

      // set the rhs values (negative vorticity)
      for (i=0;i<nx;i++) {
         for (j=0;j<ny;j++) {
            rhs[i][j] = -1.0*vort[i][j];
         }
      }

      // initialize psi
      for (i=0;i<nx;i++) {
         for (j=0;j<ny;j++) {
            psi[i][j] = 0.;
         }
      }
           //for (j=0;j<ny;j++) {
           //  psi[0][j] = 0.;
           //  psi[nx-1][j] = 0.;
           //}
           //for (i=0;i<nx;i++) {
           //  psi[i][0] = 0.;
           //  psi[i][ny-1] = 0.;
           //}

      // set BCs on psi
      if (xbdry == OPEN && ybdry == OPEN) {
        // call subroutine to solve for actual boundary velocities
        fprintf(stderr,"    solving bcs\n"); fflush(stderr);
        find_open_boundary_psi (nx,ny,vort,yf,freestream,psi);
      } else if (xbdry == OPEN) {
      //if (xbdry == OPEN) {
        // NOTE: must correct this if freestream is anything other than
        //   (1,0) or if yf != 1
        for (j=0;j<ny;j++) {
          psi[0][j] = (float)j/(float)(nx-1);
          psi[nx-1][j] = (float)j/(float)(nx-1);
        }
      } else if (ybdry == OPEN) {
        for (i=0;i<nx;i++) {
          psi[i][0] = 0.;
          psi[i][ny-1] = 1.;
        }
      }

      if (!silent) fprintf(stderr,"    running mud2sp\n"); fflush(stderr);
      mud2sp_(iparm, fparm, work, cofx, cofy,
              bndyc, rhs[0], psi[0], mgopt, &ierr);
      if (!silent) fprintf(stdout,"    mud2sp solved in %d cycles\n",iparm[16]);

      // debug print the resulting streamfunctions
      //for (i=0; i<8; i++) {
      //   for (j=0; j<8; j++) {
      //      fprintf(stdout," %g",psi[i][j]);
      //   }
      //   fprintf(stdout,"\n");
      //}

      // catch a runtime error
      if (ierr != 0) {
         fprintf(stderr,"ERROR (mud2sp_): ierr = %d\n",ierr);
         fprintf(stderr,"Quitting.\n");
         exit(0);
      }

      //exit(0);

   // Use fishpak FFT solver ----------------------------------------------
   } else {

      // set the rhs values (negative vorticity)
      for (i=0;i<nx;i++) {
         for (j=0;j<ny;j++) {
            psi[i][j] = -1.0*vort[i][j];
         }
      }

      // set the boundary conditions
      //  0 : the solution is periodic
      //  1 : the solution is specified on the boundaries
      //  3 : the derivative of the solution is specified at the boundaries
      if (xbdry == PERIODIC) {
         bdcx = 0;
      } else {
         if (isStam) {
            bdcx = 1;
         } else {
            bdcx = 1;
         }
      }
      if (ybdry == PERIODIC) {
         bdcy = 0;
      } else {
         if (isStam) {
            bdcy = 1;
         } else {
            bdcy = 1;
         }
      }

      // set boundary values (all zeroes for now)
      for (i=0;i<nx;i++) {
         bdvys[i] = 0.0;
         bdvyf[i] = 0.0;
         // vort[i][0] = 0.0;
         // vort[i][nym1] = 0.0;
      }
      //if (ybdry == WALL && !isStam) {
      if (ybdry == WALL) {
         for (i=0;i<nx;i++) {
            psi[i][0] = 0.0;
            psi[i][nym1] = 0.0;
         }
      }
      for (i=0;i<ny;i++) {
         bdvxs[i] = 0.0;
         bdvxf[i] = 0.0;
         // vort[0][i] = 0.0;
         // vort[nxm1][i] = 0.0;
      }
      //if (xbdry == WALL && !isStam) {
      if (xbdry == WALL) {
         for (i=0;i<ny;i++) {
            psi[0][i] = 0.0;
            psi[nxm1][i] = 0.0;
         }
      }

      // make the call to the external solver, note that x and y are backwards
      hwscrt_(&ys, &yf, &nym1, &bdcy, bdvys, bdvyf,
              &xs, &xf, &nxm1, &bdcx, bdvxs, bdvxf,
              &elmbda, psi[0], &nyp, &pert, &ierr, work);
      // hwscrt_(&xs, &xf, &nxm1, &bdcx, bdvxs, bdvxf,
      //         &ys, &yf, &nym1, &bdcy, bdvys, bdvyf,
      //         &elmbda, psi[0], &nxp, &pert, &ierr, work);

      // fprintf(stdout,"  hwscrt used %d locations in work array, err=%d\n",(int)(work[0]),ierr);

   // Done with if use_multigrid ------------------------------------------
   }

   // write some debug stuff
   //sprintf(outfileroot,"psi_%05d",step);
   //write_png (outfileroot,nx,ny,FALSE,FALSE,
   //                    psi,-0.1,0.2,
   //                    NULL,0.0,1.0,
   //                    NULL,0.0,1.0);

   // now, find the derivatives of streamfunction to make the velocities!
   find_gradient_of_scalar_2nd_2d (nx,ny,xbdry,ybdry,psi,NULL,v,-1.0,u,1.0);

   // if there's a mask, then find the counter-vorticity
   if (use_MASK) {

      // first, find counter-velocity
      for (i=0;i<nx;i++) {
         for (j=0;j<ny;j++) {
            //cu[i][j] = -1.0*u[i][j]*mask[i][j];
            //cv[i][j] = -1.0*v[i][j]*mask[i][j];
            cu[i][j] = -1.0*u[i][j]*(1.-mask[i][j]);
            cv[i][j] = -1.0*v[i][j]*(1.-mask[i][j]);
         }
      }

      // test-print counter-velocity
      //sprintf(outfileroot,"cu_%05d",istep);
      //write_png (outfileroot,nx,ny,FALSE,FALSE,cu,-1.0,2.0,
      //           NULL,0.0,1.0,NULL,0.0,1.0);
      //sprintf(outfileroot,"cv_%05d",istep);
      //write_png (outfileroot,nx,ny,FALSE,FALSE,cv,-1.0,2.0,
      //           NULL,0.0,1.0,NULL,0.0,1.0);

      // compute counter-vorticity
      find_curl_of_vel_2d(nx,ny,xbdry,ybdry,cu,cv,cw);

      // test-print counter-vorticity
      //sprintf(outfileroot,"cw_%05d",istep);
      //write_png (outfileroot,nx,ny,FALSE,FALSE,cw,-10.0,20.0,
      //           NULL,0.0,1.0,NULL,0.0,1.0);

      // add counter-vorticity to vorticity field and iterate!
      for (i=0;i<nx;i++) {
         for (j=0;j<ny;j++) {
            // no good: sharp points get too strong
            //vort[i][j] = vort[i][j] + cw[i][j]*(1.-mask[i][j]);
            // also no good: flat walls grow stuff...
            //vort[i][j] = vort[i][j] + cw[i][j]*mask[i][j];
            // works well, if all vort with mask>.999 are set to zero
            //vort[i][j] = vort[i][j] + cw[i][j];
            // this should be the bomb.
            //vort[i][j] = (1.-mask[i][j]) * (vort[i][j] + cw[i][j]);
            vort[i][j] = mask[i][j] * (vort[i][j] + cw[i][j]);
         }
      }

      // finally, zero out vorticity fully inside of mask
      //for (i=0;i<nx;i++) {
      //   for (j=0;j<ny;j++) {
      //      if (mask[i][j] > 0.999) vort[i][j] = 0.;
      //   }
      //}

      // check difference between this and last iterations - completion check

      if (istep > 0) {
         // first, find maximum vorticity
         maxvort = 0.;
         for (i=0;i<nx;i++) {
            for (j=0;j<ny;j++) {
               temp = fabs(vort[i][j]);
               if (temp > maxvort) maxvort = temp;
            }
         }
         // now, find maximum difference between this and previous counter vorticity
         maxcorr = 0.;
         for (i=0;i<nx;i++) {
            for (j=0;j<ny;j++) {
               temp = fabs(cw_old[i][j] - cw[i][j]);
               if (temp > maxcorr) maxcorr = temp;
            }
         }
         // if the error is small, jump out of the loop
         fprintf(stderr,"  it %d, err %g\n",istep,maxcorr/maxvort);
         if (maxcorr/maxvort < maskerr) {
            // yay! we're done iterating
            istep += maxstep;
            // reset back to full solve
            fparm[4] = muderr;
            iparm[12] = 100;
         } else if (maxcorr/maxvort < 5.*maskerr) {
            // always run a specific number of cycles
            fparm[4] = 0.0;
            iparm[12] = 2;
         } else {
            // always run a specific number of cycles
            fparm[4] = 0.0;
            iparm[12] = 1;
         }
         if (isnan(maxcorr/maxvort)) {
            if (maxcorr < 1.e-12) {
               // yay! we're done iterating
               istep += maxstep;
               // reset back to full solve
               fparm[4] = muderr;
               iparm[12] = 100;
            } else {
               // whoah, something bad happened
               fprintf(stderr,"Quitting.\n");
               exit(1);
            }
         }
      }
      // copy over the current to the old
      for (i=0;i<nx;i++) {
         for (j=0;j<ny;j++) {
            cw_old[i][j] = cw[i][j];
         }
      }

      //fprintf(stderr,"  iteration step %d\n",istep);
   }

   // try forcing the driven cavity this way!
   if (driven_cavity) {
      for (i=0;i<nx;i++) {
         u[i][nym1] = 1.0;
      }
   }


   }  // end iterative solution
   //exit(0);

   return(0);
}


/*
 * find_open_boundary_psi takes a vorticity field, and solves for the velocity field
         xbdry = OPEN;
         ybdry = OPEN;
 *
 */
int find_open_boundary_psi (int nx,int ny,float **vort,
      float yf,float *freestream,float **psi) {

   int i,j;
   float xs = 0.0;
   float ys = 0.0;
   float xf = 1.0;
   //float yf = SCALE;
   float dx = 1./(nx-1);
   float dy = yf/(ny-1);
   float thisvel[2];
   float *vel;

   // malloc space for velocities
   if (nx>ny) vel = allocate_1d_array_f((long int)nx);
   else vel = allocate_1d_array_f((long int)ny);

   // initialize the static variables
   find_biot_savart (TRUE,0.,0.,nx,ny,dx,dy,vort,thisvel);

   // first, march across bottom row (j=0, y=0)
   #pragma omp parallel for private(i,thisvel)
   for (i=0;i<nx;i++) {
     find_biot_savart (FALSE,i*dx,0.,nx,ny,dx,dy,vort,thisvel);
     vel[i] = freestream[1] + thisvel[1];
   }
   //exit(0);
   psi[0][0] = 0.0;
   for (i=1;i<nx;i++) {
     psi[i][0] = psi[i-1][0] - 0.5*(vel[i]+vel[i-1])*dx;
     //fprintf(stderr,"x=%g  v=%g  psi=%g\n",i*dx,vel[i],psi[i][0]);
   }

   // then march up left side
   #pragma omp parallel for private(j,thisvel)
   for (j=0;j<ny;j++) {
     find_biot_savart (FALSE,0.,j*dy,nx,ny,dx,dy,vort,thisvel);
     vel[j] = freestream[0] + thisvel[0];
   }
   for (j=1;j<ny;j++) {
     psi[0][j] = psi[0][j-1] + 0.5*(vel[j]+vel[j-1])*dy;
     //fprintf(stderr,"y=%g  u=%g  psi=%g\n",j*dy,vel[j],psi[0][j]);
   }

   // then march up right side
   #pragma omp parallel for private(j,thisvel)
   for (j=0;j<ny;j++) {
     find_biot_savart (FALSE,1.,j*dy,nx,ny,dx,dy,vort,thisvel);
     vel[j] = freestream[0] + thisvel[0];
   }
   for (j=1;j<ny;j++) {
     psi[nx-1][j] = psi[nx-1][j-1] + 0.5*(vel[j]+vel[j-1])*dy;
     //fprintf(stderr,"y=%g  u=%g  psi=%g\n",j*dy,vel[j],psi[nx-1][j]);
   }

   // finally, march across top row (j=ny-1, y=yf)
   #pragma omp parallel for private(i,thisvel)
   for (i=0;i<nx;i++) {
     find_biot_savart (FALSE,i*dx,yf,nx,ny,dx,dy,vort,thisvel);
     vel[i] = freestream[1] + thisvel[1];
   }
   for (i=1;i<nx;i++) {
     psi[i][ny-1] = psi[i-1][ny-1] - 0.5*(vel[i]+vel[i-1])*dx;
     //fprintf(stderr,"x=%g  v=%g  psi=%g\n",i*dx,vel[i],psi[i][ny-1]);
   }
   //exit(0);
   free_1d_array_f(vel);

   return(0);
}

int find_biot_savart (const int initialize,
                      const float xp,const float yp,
                      const int nx,const int ny,
                      const float dx,const float dy,
                      float **vort,float *vel) {

  const int treecode = FALSE;
  static int istart = 100000;
  static int iend = -1;
  static int jstart = 100000;
  static int jend = -1;
  float xdist,ydist,distsq;
  static int nlevels = -1;
  static float** va[20];
  static int firsttime = TRUE;

  // initialize some arrays
  if (initialize) {

    istart = 100000;
    iend = -1;
    jstart = 100000;
    jend = -1;
    
    // which columns have non-zero circulation?
    // which rows have non-zero circulation?
    for (int i=0;i<nx;i++) {
      for (int j=0;j<ny;j++) {
        if (abs(vort[i][j]) > 1.e-20) {
          if (i<istart) istart=i;
          if (i+1>iend) iend=i+1;
          if (j<jstart) jstart=j;
          if (j+1>jend) jend=j+1;
        }
      }
    }

    // for O(NlogN) version: generate the hierarchical array
    if (treecode) {
      int ilevel = 0;
      va[ilevel] = vort;
      int oldnx = nx;
      int oldny = ny;
      while (oldnx > 4 && oldny > 4 && ilevel < 20) {
        ilevel++;
        int newnx = (oldnx+1)/2;
        int newny = (oldny+1)/2;
        if (firsttime) va[ilevel] = allocate_2d_array_f(newnx, newny);
        for (int i=0; i<newnx; i++) for (int j=0; j<newny; j++) va[ilevel][i][j] = 0.;
        for (int i=0; i<oldnx; i++) {
          for (int j=0; j<oldny; j++) {
            va[ilevel][i/2][j/2] += va[ilevel-1][i][j];
          }
        }
        iend = newnx;
        jend = newny;
        oldnx = newnx;
        oldny = newny;
      }
      nlevels = ilevel;
    }

    // never alloc again
    if (firsttime) firsttime = FALSE;
    //fprintf(stderr,"%d %d  %d %d\n",istart,iend,jstart,jend);
  }

  // now, only iterate over the necessary rows/columns
  vel[0] = 0.;
  vel[1] = 0.;
  if (treecode) {
    // iterate over coarsest grid first
    for (int i=0;i<iend;i++) {
      for (int j=0;j<jend;j++) {
        recursive_biot_savart(xp,yp,i,j,nlevels,
           dx*pow(2,nlevels),dy*pow(2,nlevels),
           va,vel);
      }
    }

  } else {  // direct method
    for (int i=istart;i<iend;i++) {
      xdist = i*dx - xp;
      for (int j=jstart;j<jend;j++) {
        if (abs(vort[i][j]) > 1.e-20) {
          ydist = j*dy - yp;
          distsq = xdist*xdist+ydist*ydist+dx*dy;
          vel[0] += ydist * vort[i][j] / distsq;
          vel[1] -= xdist * vort[i][j] / distsq;
        }
      }
    }
  }
  vel[0] *= dx*dy/(2.*M_PI);
  vel[1] *= dx*dy/(2.*M_PI);

  return ((iend-istart-1)*(jend-jstart-1));
}


/*
 * recursive Biot-Savart law, for hierarchical treecode
 *
 * needs work: fix distances to account for offset box centers
 *             prevent recursing into oob i,j
 */
void recursive_biot_savart (const float xp, const float yp,
      const int i, const int j, const int level,
      const float dx, const float dy,
      float*** va, float* vel) {

   // don't even bother if vorticity is weak
   if (abs(va[level][i][j]) < 1.e-20) return;

   // check box-opening criterion
   const float xdist = i*dx - xp;
   const float ydist = j*dy - yp;
   float distsq = xdist*xdist+ydist*ydist;
   if (level == 0 || distsq > 3.*dx) {
      // compute influence at this level
      vel[0] += ydist * va[level][i][j] / distsq;
      vel[1] -= xdist * va[level][i][j] / distsq;

   } else {
      // recurse to finer level
      recursive_biot_savart(xp,yp,2*i,2*j,level-1,dx/2.,dy/2.,va,vel);
      recursive_biot_savart(xp,yp,2*i,2*j+1,level-1,dx/2.,dy/2.,va,vel);
      recursive_biot_savart(xp,yp,2*i+1,2*j,level-1,dx/2.,dy/2.,va,vel);
      recursive_biot_savart(xp,yp,2*i+1,2*j+1,level-1,dx/2.,dy/2.,va,vel);
   }

   return;
}

/*
 * differentiate the a scalar to find the velocity
 *
 * NOT USED ANY MORE, see find_gradient_of_scalar_2nd_2d
 */
int find_gradient_of_scalar_2d (int nx,int ny,int xbdry,int ybdry,float **psi,float **u,float umult,float **v,float vmult) {

   int i,j;
   int nxm1 = nx-1;
   int nym1 = ny-1;
   int im1,ip1,jm1,jp1;
   int i0im1,inip1;
   int j0jm1,jnjp1;
   FLOAT i0cx,j0cy;
   FLOAT incx,jncy;
   FLOAT hxi,hyi,cx,cy;
   int ini,jnj,ii,jj;
   hxi = 0.5*umult*(nxm1)/SCALE;		// negative because vel is neg gradient
   // hyi = 0.5*(nym1);
   hyi = 0.5*vmult*(nxm1);

   // now, account for periodic BCs by setting some indexes
   if (xbdry == WALL) {
      i0im1 = 0;
      i0cx = 2.0;
      ini = nxm1;
      inip1 = nxm1;
      incx = 2.0;
   } else if (xbdry == PERIODIC) {
      i0im1 = nxm1 - 1;
      i0cx = 1.0;
      ini = 0;
      inip1 = 1;
      incx = 1.0;
   }
 
   if (ybdry == WALL) {
      j0jm1 = 0;
      j0cy = 2.0;
      jnj = nym1;
      jnjp1 = nym1;
      jncy = 2.0;
   } else if (ybdry == PERIODIC) {
      j0jm1 = nym1 - 1;
      j0cy = 1.0;
      jnj = 0;
      jnjp1 = 1;
      jncy = 1.0;
   }

   for (i=0; i<nx; i++) {
      im1 = i-1;
      ii = i;
      ip1 = i+1;
      cx = 1.0;
      if (i == 0) {
         im1 = i0im1;
         ii = i;
         ip1 = 1;
         cx = i0cx;
      } else if (i == nxm1-1) {
         // im1 is correct
         ip1 = ini;
      } else if (i == nxm1) {
         im1 = nxm1 - 1;
         ii = ini;
         ip1 = inip1;
         cx = incx;
      }
 
   for (j=0; j<ny; j++) {
      jm1 = j-1;
      jj = j;
      jp1 = j+1;
      cy = 1.0;
      if (j == 0) {
         jm1 = j0jm1;
         jp1 = 1;
         cy = j0cy;
      } else if (j == nym1-1) {
         // jm1 is correct
         jp1 = jnj;
      } else if (j == nym1) {
         jm1 = nym1 - 1;
         jj = jnj;
         jp1 = jnjp1;
         cy = jncy;
      }

      // partial in x-direction
      u[i][j] = (psi[ip1][jj]-psi[im1][jj])*hxi*cx;
      // partial in y-direction
      v[i][j] = (psi[ii][jp1]-psi[ii][jm1])*hyi*cy;
 
   }
   }
   // fprintf(stdout,"  vel at center is %g %g\n",u[50][50],v[50][50]);

   return(0);
}


/*
 * differentiate the a scalar to find the velocity
 *
 * use a new method that is 2nd order everywhere
 * and if mask is not NULL, zeros gradients where mask is present
 */
int find_gradient_of_scalar_2nd_2d (int nx,int ny,int xbdry,int ybdry,float **psi,float **mask,float **u,float umult,float **v,float vmult) {

   int i,j;
   int nxm1 = nx-1;
   int nym1 = ny-1;
   FLOAT hxi,hyi;

   // new school, max domain is 0..1
   if (nx > ny) {
      hxi = 0.5*umult*(nxm1);
      hyi = 0.5*vmult*(nxm1);
   } else {
      hxi = 0.5*umult*(nym1);
      hyi = 0.5*vmult*(nym1);
   }

   // do the middle first
   #pragma omp parallel for private(i,j)
   for (i=1; i<nxm1; i++) {
   for (j=1; j<nym1; j++) {
      // partial in x-direction
      u[i][j] = (psi[i+1][j]-psi[i-1][j])*hxi;
      // partial in y-direction
      v[i][j] = (psi[i][j+1]-psi[i][j-1])*hyi;
   }
   }

   // now, do left and right sides
   if (xbdry == WALL || xbdry == OPEN) {
      i = 0;
      for (j=1; j<nym1; j++) {
         u[i][j] = (-1.*psi[i+2][j]+4.*psi[i+1][j]-3.*psi[i][j])*hxi;
         v[i][j] = (psi[i][j+1]-psi[i][j-1])*hyi;
      }
      i = nxm1;
      for (j=1; j<nym1; j++) {
         u[i][j] = (psi[i-2][j]-4.*psi[i-1][j]+3.*psi[i][j])*hxi;
         v[i][j] = (psi[i][j+1]-psi[i][j-1])*hyi;
      }
   } else if (xbdry == PERIODIC) {
      i = 0;
      for (j=1; j<nym1; j++) {
         u[i][j] = (psi[i+1][j]-psi[nxm1-1][j])*hxi;
         v[i][j] = (psi[i][j+1]-psi[i][j-1])*hyi;
         u[nxm1][j] = u[i][j];
         v[nxm1][j] = v[i][j];
      }
   }
 
   // now, do top and bottom sides
   if (ybdry == WALL || ybdry == OPEN) {
      j = 0;
      for (i=1; i<nxm1; i++) {
         u[i][j] = (psi[i+1][j]-psi[i-1][j])*hxi;
         v[i][j] = (-1.*psi[i][j+2]+4.*psi[i][j+1]-3.*psi[i][j])*hyi;
      }
      j = nym1;
      for (i=1; i<nxm1; i++) {
         u[i][j] = (psi[i+1][j]-psi[i-1][j])*hxi;
         v[i][j] = (psi[i][j-2]-4.*psi[i][j-1]+3.*psi[i][j])*hyi;
      }
   } else if (ybdry == PERIODIC) {
      j = 0;
      for (i=1; i<nxm1; i++) {
         u[i][j] = (psi[i+1][j]-psi[i-1][j])*hxi;
         v[i][j] = (psi[i][j+1]-psi[i][nym1-1])*hyi;
         u[i][nym1] = u[i][j];
         v[i][nym1] = v[i][j];
      }
   }

   // finally, do the corners
   if ((xbdry == WALL && ybdry == WALL) || (xbdry == OPEN && ybdry == OPEN)) {
      i = 0; j = 0;
      u[i][j] = (-1.*psi[i+2][j]+4.*psi[i+1][j]-3.*psi[i][j])*hxi;
      v[i][j] = (-1.*psi[i][j+2]+4.*psi[i][j+1]-3.*psi[i][j])*hyi;
      i = 0; j = nym1;
      u[i][j] = (-1.*psi[i+2][j]+4.*psi[i+1][j]-3.*psi[i][j])*hxi;
      v[i][j] = (psi[i][j-2]-4.*psi[i][j-1]+3.*psi[i][j])*hyi;
      i = nxm1; j = nym1;
      u[i][j] = (psi[i-2][j]-4.*psi[i-1][j]+3.*psi[i][j])*hxi;
      v[i][j] = (psi[i][j-2]-4.*psi[i][j-1]+3.*psi[i][j])*hyi;
      i = nxm1; j = 0;
      u[i][j] = (psi[i-2][j]-4.*psi[i-1][j]+3.*psi[i][j])*hxi;
      v[i][j] = (-1.*psi[i][j+2]+4.*psi[i][j+1]-3.*psi[i][j])*hyi;
   } else if (xbdry == PERIODIC && ybdry == WALL) {
      i = 0; j = 0;
      u[i][j] = (psi[i+1][j]-psi[nxm1-1][j])*hxi;
      v[i][j] = (-1.*psi[i][j+2]+4.*psi[i][j+1]-3.*psi[i][j])*hyi;
      i = 0; j = nym1;
      u[i][j] = (psi[i+1][j]-psi[nxm1-1][j])*hxi;
      v[i][j] = (psi[i][j-2]-4.*psi[i][j-1]+3.*psi[i][j])*hyi;
      i = nxm1; j = nym1;
      u[i][j] = u[0][j];
      v[i][j] = v[0][j];
      i = nxm1; j = 0;
      u[i][j] = u[0][j];
      v[i][j] = v[0][j];
   } else if (xbdry == WALL && ybdry == PERIODIC) {
      i = 0; j = 0;
      u[i][j] = (-1.*psi[i+2][j]+4.*psi[i+1][j]-3.*psi[i][j])*hxi;
      v[i][j] = (psi[i][j+1]-psi[i][nym1-1])*hyi;
      i = nxm1; j = 0;
      u[i][j] = (psi[i-2][j]-4.*psi[i-1][j]+3.*psi[i][j])*hxi;
      v[i][j] = (psi[i][j+1]-psi[i][nym1-1])*hyi;
      i = 0; j = nym1;
      u[i][j] = u[i][0];
      v[i][j] = v[i][0];
      i = nxm1; j = nym1;
      u[i][j] = u[i][0];
      v[i][j] = v[i][0];
   } else {
      i = 0; j = 0;
      u[i][j] = (psi[i+1][j]-psi[nxm1-1][j])*hxi;
      v[i][j] = (psi[i][j+1]-psi[i][nym1-1])*hyi;
      i = nxm1; j = 0;
      u[i][j] = u[0][0];
      v[i][j] = v[0][0];
      i = 0; j = nym1;
      u[i][j] = u[0][0];
      v[i][j] = v[0][0];
      i = nxm1; j = nym1;
      u[i][j] = u[0][0];
      v[i][j] = v[0][0];
   }

   // if there is a mask, many gradients will become zero
   if (mask != NULL) {

      // first, the cell itself
      #pragma omp parallel for private(i,j)
      for (i=0; i<nx; i++) {
      for (j=0; j<ny; j++) {
         // partial in x-direction
         u[i][j] *= mask[i][j];
         // partial in y-direction
         v[i][j] *= mask[i][j];
      }
      }

      // the cell left
      #pragma omp parallel for private(i,j)
      for (i=1; i<nx; i++) {
      for (j=0; j<ny; j++) {
         // partial in x-direction
         u[i][j] *= mask[i-1][j];
      }
      }

      // the cell right
      #pragma omp parallel for private(i,j)
      for (i=0; i<nxm1; i++) {
      for (j=0; j<ny; j++) {
         // partial in x-direction
         u[i][j] *= mask[i+1][j];
      }
      }

      // the cell above
      #pragma omp parallel for private(i,j)
      for (i=0; i<nx; i++) {
      for (j=0; j<nym1; j++) {
         // partial in y-direction
         v[i][j] *= mask[i][j+1];
      }
      }

      // the cell below
      #pragma omp parallel for private(i,j)
      for (i=0; i<nx; i++) {
      for (j=1; j<ny; j++) {
         // partial in y-direction
         v[i][j] *= mask[i][j-1];
      }
      }

   }

   return(0);
}


/*
 * calculate the curl of the velocity field
 */
int find_curl_of_vel_2d(int nx,int ny,int xbdry,int ybdry,
      float **u,float **v,float **vort) {

   int i,j;
   int nxm1 = nx-1;
   int nym1 = ny-1;
   FLOAT hxi,hyi;

   // old way
   hxi = 0.5*(nxm1)/SCALE;
   hyi = 0.5*(nxm1);		// and *nxm1 because xsize=1
   // new way
   if (nx > ny) {
      hxi = 0.5*(nxm1);
      hyi = 0.5*(nxm1);
   } else {
      hxi = 0.5*(nym1);
      hyi = 0.5*(nym1);
   }

   // do the middle first
   #pragma omp parallel for private(i,j)
   for (i=1; i<nxm1; i++) {
   for (j=1; j<nym1; j++) {
      vort[i][j] = (v[i+1][j]-v[i-1][j])*hxi -
                   (u[i][j+1]-u[i][j-1])*hyi;
   }
   }

   // now, do left and right sides
   if (xbdry == WALL || xbdry == OPEN) {
      i = 0;
      for (j=1; j<nym1; j++) {
         vort[i][j] = (-1.*v[i+2][j]+4.*v[i+1][j]-3.*v[i][j])*hxi -
                      (u[i][j+1]-u[i][j-1])*hyi;
      }
      i = nxm1;
      for (j=1; j<nym1; j++) {
         vort[i][j] = (v[i-2][j]-4.*v[i-1][j]+3.*v[i][j])*hxi -
                      (u[i][j+1]-u[i][j-1])*hyi;
      }
   } else if (xbdry == PERIODIC) {
      i = 0;
      for (j=1; j<nym1; j++) {
         vort[i][j] = (v[i+1][j]-v[nxm1-1][j])*hxi -
                      (u[i][j+1]-u[i][j-1])*hyi;
         vort[nxm1][j] = vort[i][j];
      }
   }
 
   // now, do top and bottom sides
   if (ybdry == WALL || ybdry == OPEN) {
      j = 0;
      for (i=1; i<nxm1; i++) {
         vort[i][j] = (v[i+1][j]-v[i-1][j])*hxi - 
                      (-1.*u[i][j+2]+4.*u[i][j+1]-3.*u[i][j])*hyi;
      }
      j = nym1;
      for (i=1; i<nxm1; i++) {
         vort[i][j] = (v[i+1][j]-v[i-1][j])*hxi - 
                      (u[i][j-2]-4.*u[i][j-1]+3.*u[i][j])*hyi;
      }
   } else if (ybdry == PERIODIC) {
      j = 0;
      for (i=1; i<nxm1; i++) {
         vort[i][j] = (v[i+1][j]-v[i-1][j])*hxi -
                      (u[i][j+1]-u[i][nym1-1])*hyi;
         vort[i][nym1] = vort[i][j];
      }
   }

   // finally, do the corners
   if ((xbdry == WALL && ybdry == WALL) || (xbdry == OPEN && ybdry == OPEN)) {
      i = 0; j = 0;
      vort[i][j] = (-1.*v[i+2][j]+4.*v[i+1][j]-3.*v[i][j])*hxi -
                   (-1.*u[i][j+2]+4.*u[i][j+1]-3.*u[i][j])*hyi;
      i = 0; j = nym1;
      vort[i][j] = (-1.*v[i+2][j]+4.*v[i+1][j]-3.*v[i][j])*hxi -
                   (u[i][j-2]-4.*u[i][j-1]+3.*u[i][j])*hyi;
      i = nxm1; j = nym1;
      vort[i][j] = (v[i-2][j]-4.*v[i-1][j]+3.*v[i][j])*hxi -
                   (u[i][j-2]-4.*u[i][j-1]+3.*u[i][j])*hyi;
      i = nxm1; j = 0;
      vort[i][j] = (v[i-2][j]-4.*v[i-1][j]+3.*v[i][j])*hxi -
                   (-1.*u[i][j+2]+4.*u[i][j+1]-3.*u[i][j])*hyi;
   } else if (xbdry == PERIODIC && ybdry == WALL) {
      i = 0; j = 0;
      vort[i][j] = (v[i+1][j]-v[nxm1-1][j])*hxi -
                   (-1.*u[i][j+2]+4.*u[i][j+1]-3.*u[i][j])*hyi;
      i = 0; j = nym1;
      vort[i][j] = (v[i+1][j]-v[nxm1-1][j])*hxi -
                   (u[i][j-2]-4.*u[i][j-1]+3.*u[i][j])*hyi;
      i = nxm1; j = nym1;
      vort[i][j] = vort[0][j];
      i = nxm1; j = 0;
      vort[i][j] = vort[0][j];
   } else if (xbdry == WALL && ybdry == PERIODIC) {
      i = 0; j = 0;
      vort[i][j] = (-1.*v[i+2][j]+4.*v[i+1][j]-3.*v[i][j])*hxi -
                   (u[i][j+1]-u[i][nym1-1])*hyi;
      i = nxm1; j = 0;
      vort[i][j] = (v[i-2][j]-4.*v[i-1][j]+3.*v[i][j])*hxi -
                   (u[i][j+1]-u[i][nym1-1])*hyi;
      i = 0; j = nym1;
      vort[i][j] = vort[i][0];
      i = nxm1; j = nym1;
      vort[i][j] = vort[i][0];
   } else { // xbdry == PERIODIC && ybdry == PERIODIC
      i = 0; j = 0;
      vort[i][j] = (v[i+1][j]-v[nxm1-1][j])*hxi -
                   (u[i][j+1]-u[i][nym1-1])*hyi;
      i = nxm1; j = 0;
      vort[i][j] = vort[0][0];
      i = 0; j = nym1;
      vort[i][j] = vort[0][0];
      i = nxm1; j = nym1;
      vort[i][j] = vort[0][0];
   }

   return(0);
}


/*
 * moc_advect_2d is the implicit scheme for updating all field values
 */
int moc_advect_2d (int nx,int ny,int xbdry,int ybdry,
      float** mask,float** u,float** v,float* fs,
      float*** in,float*** out,
      float dt,int order,int moveColors) {

   enum interpMeth {
      cic, tsc, m4p, cic2, cic3 } interp = m4p;
   int i,j,k;
   int sc_cnt;
   float xf,yf,dx,accx,accy,oodt;
   float px,py,newx,newy,temp;
   float u0,u1,u2,u3,v0,v1,v2,v3;
   float outvals[MAX_SCALARS];
   float **tempin[MAX_SCALARS];
   float **tempout[MAX_SCALARS];

   oodt = 1.0/dt;

   // create temporary array pointing to all non-NULL in/out arrays!
   // this prevents the need for many checks vs. NULL
   sc_cnt = 0;
   for (i=0; i<MAX_SCALARS; i++) {
      tempin[i] = NULL;
      tempout[i] = NULL;
      if (in[i] != NULL) {
         tempin[sc_cnt] = in[i];
         tempout[sc_cnt] = out[i];
         sc_cnt++;
      }
   }

   // fprintf(stderr,"in moc_advect_2d\n"); fflush(stderr);
   if (nx > ny) {
      xf = 1.0;
      yf = (float)(ny-1)/(float)(nx-1);
      dx = 1.0 / (float)(nx-1);
   } else {
      xf = (float)(nx-1)/(float)(ny-1);
      yf = 1.0;
      dx = 1.0 / (float)(ny-1);
   }

   // for each point, march backwards in space and time to find the value
   if (order == 1) {
      // use a 1-step method
      #pragma omp parallel for private(i,px,j,py,newx,newy,k,outvals)
      for (i=0;i<nx;i++) {
         px = (float)i * dx;
         for (j=0;j<ny;j++) {
            int do_this_pixel = TRUE;
            float velmult = 1.;
            if (mask != NULL) {
               if (mask[i][j] < 1.e-5) do_this_pixel = FALSE;
               else velmult = mask[i][j];
            }
            if (i==0 && fs[0] > 0.0) do_this_pixel = FALSE;
            if (i==nx-1 && fs[0] < 0.0) do_this_pixel = FALSE;
            if (j==0 && fs[1] > 0.0) do_this_pixel = FALSE;
            if (j==ny-1 && fs[1] < 0.0) do_this_pixel = FALSE;

            if (do_this_pixel) {
               py = (float)j * dx;
               newx = px-dt*velmult*u[i][j];
               newy = py-dt*velmult*v[i][j];
               //if (xbdry == WALL || xbdry == OPEN) {
               //   if (newx > xf-EPSILON) newx = xf-EPSILON;
               //   if (newx < 0.0+EPSILON) newx = EPSILON;
               //}
               //if (ybdry == WALL || ybdry == OPEN) {
               //   if (newy > yf-EPSILON) newy = yf-EPSILON;
               //   if (newy < 0.0+EPSILON) newy = EPSILON;
               //}
               if (interp == cic)
                  interpolate_array_using_CIC_2d(nx,ny,xbdry,ybdry,tempin,newx,newy,sc_cnt,outvals);
               else if (interp == tsc)
                  interpolate_array_using_TSC_2d(nx,ny,xbdry,ybdry,tempin,newx,newy,sc_cnt,outvals);
               else
                  interpolate_array_using_M4p_2d(nx,ny,xbdry,ybdry,mask,tempin,newx,newy,sc_cnt,outvals);
               for (k=0; k<sc_cnt; k++) tempout[k][i][j] = outvals[k];
            } else {
               for (k=0; k<sc_cnt; k++) tempout[k][i][j] = tempin[k][i][j];
            }
         }
      }
   } else if (order == 2) {
      // use a two-step method
      #pragma omp parallel for private(i,px,j,py,u0,u1,v0,v1,newx,newy,k,outvals)
      for (i=0;i<nx;i++) {
         px = (float)i * dx;
         for (j=0;j<ny;j++) {
            int do_this_pixel = TRUE;
            float velmult = 1.;
            if (mask != NULL) {
               if (mask[i][j] < 1.e-5) do_this_pixel = FALSE;
               else velmult = mask[i][j];
            }
            if (i==0 && fs[0] > 0.0) do_this_pixel = FALSE;
            if (i==nx-1 && fs[0] < 0.0) do_this_pixel = FALSE;
            if (j==0 && fs[1] > 0.0) do_this_pixel = FALSE;
            if (j==ny-1 && fs[1] < 0.0) do_this_pixel = FALSE;

            if (do_this_pixel) {
               py = (float)j * dx;
               // find vel at px,py
               u0 = u[i][j];
               v0 = v[i][j];

               // find position 1 explicit Euler step backwards
               newx = px-dt*velmult*u0;
               //if (xbdry == WALL || xbdry == OPEN) {
               //   if (newx > xf-EPSILON) newx = xf-EPSILON;
               //   if (newx < 0.0+EPSILON) newx = EPSILON;
               //}
               newy = py-dt*velmult*v0;
               //if (ybdry == WALL || ybdry == OPEN) {
               //   if (newy > yf-EPSILON) newy = yf-EPSILON;
               //   if (newy < 0.0+EPSILON) newy = EPSILON;
               //}

               // find vel at newx,newy
               if (interp == cic)
                  interpolate_vel_using_CIC_2d(nx,ny,xbdry,ybdry,u,v,newx,newy,&u1,&v1);
               else if (interp == tsc)
                  interpolate_vel_using_TSC_2d(nx,ny,xbdry,ybdry,u,v,newx,newy,&u1,&v1);
               else
                  interpolate_vel_using_M4p_2d(nx,ny,xbdry,ybdry,mask,u,v,newx,newy,&u1,&v1);

               // find position back 1 step using average of two velocities
               newx = px-dt*0.5*velmult*(u0+u1);
               //if (xbdry == WALL || xbdry == OPEN) {
               //   if (newx > xf-EPSILON) newx = xf-EPSILON;
               //   if (newx < 0.0+EPSILON) newx = EPSILON;
               //}
               newy = py-dt*0.5*velmult*(v0+v1);
               //if (ybdry == WALL || ybdry == OPEN) {
               //   if (newy > yf-EPSILON) newy = yf-EPSILON;
               //   if (newy < 0.0+EPSILON) newy = EPSILON;
               //}

               // are we able to calculate the acceleration here?

               // finally, return scalar-valued vorticity from there
               if (interp == cic)
                  interpolate_array_using_CIC_2d(nx,ny,xbdry,ybdry,tempin,newx,newy,sc_cnt,outvals);
               else if (interp == tsc)
                  interpolate_array_using_TSC_2d(nx,ny,xbdry,ybdry,tempin,newx,newy,sc_cnt,outvals);
               else
                  interpolate_array_using_M4p_2d(nx,ny,xbdry,ybdry,mask,tempin,newx,newy,sc_cnt,outvals);
               for (k=0; k<sc_cnt; k++) tempout[k][i][j] = outvals[k];
            } else {
               for (k=0; k<sc_cnt; k++) tempout[k][i][j] = tempin[k][i][j];
            }
         }
      }
   } else if (order == 4) {
      // use a four-step method, RK4
      #pragma omp parallel for private(i,px,j,py,u0,u1,u2,u3,v0,v1,v2,v3,newx,newy,accx,accy,k,outvals)
      for (i=0;i<nx;i++) {
         px = (float)i * dx;
         for (j=0;j<ny;j++) {
            int do_this_pixel = TRUE;
            float velmult = 1.;
            if (mask != NULL) {
               if (mask[i][j] < 1.e-5) do_this_pixel = FALSE;
               else velmult = mask[i][j];
            }
            if (i==0 && fs[0] > 0.0) do_this_pixel = FALSE;
            if (i==nx-1 && fs[0] < 0.0) do_this_pixel = FALSE;
            if (j==0 && fs[1] > 0.0) do_this_pixel = FALSE;
            if (j==ny-1 && fs[1] < 0.0) do_this_pixel = FALSE;

            if (do_this_pixel) {
               py = (float)j * dx;
               // find vel at px,py (k_1)
               u0 = u[i][j];
               v0 = v[i][j];

               // find position 1 explicit Euler step backwards
               newx = px-dt*0.5*velmult*u0;
               //if (xbdry == WALL || xbdry == OPEN) {
               //   if (newx > xf-EPSILON) newx = xf-EPSILON;
               //   if (newx < 0.0+EPSILON) newx = EPSILON;
               //}
               newy = py-dt*0.5*velmult*v0;
               //if (ybdry == WALL || ybdry == OPEN) {
               //   if (newy > yf-EPSILON) newy = yf-EPSILON;
               //   if (newy < 0.0+EPSILON) newy = EPSILON;
               //}
               // find vel at newx,newy (k_2)
               if (interp == cic)
                  interpolate_vel_using_CIC_2d(nx,ny,xbdry,ybdry,u,v,newx,newy,&u1,&v1);
               else if (interp == tsc)
                  interpolate_vel_using_TSC_2d(nx,ny,xbdry,ybdry,u,v,newx,newy,&u1,&v1);
               else
                  interpolate_vel_using_M4p_2d(nx,ny,xbdry,ybdry,mask,u,v,newx,newy,&u1,&v1);

               // find position 2 explicit Euler step backwards
               newx = px-dt*0.5*velmult*u1;
               //if (xbdry == WALL || xbdry == OPEN) {
               //   if (newx > xf-EPSILON) newx = xf-EPSILON;
               //   if (newx < 0.0+EPSILON) newx = EPSILON;
               //}
               newy = py-dt*0.5*velmult*v1;
               //if (ybdry == WALL || ybdry == OPEN) {
               //   if (newy > yf-EPSILON) newy = yf-EPSILON;
               //   if (newy < 0.0+EPSILON) newy = EPSILON;
               //}
               // find vel at newx,newy (k_3)
               if (interp == cic)
                  interpolate_vel_using_CIC_2d(nx,ny,xbdry,ybdry,u,v,newx,newy,&u2,&v2);
               else if (interp == tsc)
                  interpolate_vel_using_TSC_2d(nx,ny,xbdry,ybdry,u,v,newx,newy,&u2,&v2);
               else
                  interpolate_vel_using_M4p_2d(nx,ny,xbdry,ybdry,mask,u,v,newx,newy,&u2,&v2);

               // find position 3 explicit Euler step backwards
               newx = px-dt*velmult*u2;
               //if (xbdry == WALL || xbdry == OPEN) {
               //   if (newx > xf-EPSILON) newx = xf-EPSILON;
               //   if (newx < 0.0+EPSILON) newx = EPSILON;
               //}
               newy = py-dt*velmult*v2;
               //if (ybdry == WALL || ybdry == OPEN) {
               //   if (newy > yf-EPSILON) newy = yf-EPSILON;
               //   if (newy < 0.0+EPSILON) newy = EPSILON;
               //}
               // find vel at newx,newy (k_4)
               if (interp == cic)
                  interpolate_vel_using_CIC_2d(nx,ny,xbdry,ybdry,u,v,newx,newy,&u3,&v3);
               else if (interp == tsc)
                  interpolate_vel_using_TSC_2d(nx,ny,xbdry,ybdry,u,v,newx,newy,&u3,&v3);
               else
                  interpolate_vel_using_M4p_2d(nx,ny,xbdry,ybdry,mask,u,v,newx,newy,&u3,&v3);

               // find position back 1 step using average of four velocities
               accx = 0.16666667*(u0+u3+2.*(u1+u2));
               newx = px-dt*velmult*accx;
               //if (xbdry == WALL || xbdry == OPEN) {
               //   if (newx > xf-EPSILON) newx = xf-EPSILON;
               //   if (newx < 0.0+EPSILON) newx = EPSILON;
               //}
               accy = 0.16666667*(v0+v3+2.*(v1+v2));
               newy = py-dt*velmult*accy;
               //if (ybdry == WALL || ybdry == OPEN) {
               //   if (newy > yf-EPSILON) newy = yf-EPSILON;
               //   if (newy < 0.0+EPSILON) newy = EPSILON;
               //}

               // gather an estimate of the acceleration (u1-u0)/dt?
               accx = velmult*(u0-accx)*oodt;
               accy = velmult*(v0-accy)*oodt;
               //accmag = pow(accx,2) + pow(accy,2);
               //accmag = sqrt(accmag);
               //fprintf(stdout,"%d %d %g\n",i,j,accmag);
               //fprintf(stdout,"%d %d %g %g\n",i,j,accx,accy);
               // somehow, accx looks like accel in y, and vice versa

               // finally, return scalar values from there
               if (interp == cic)
                  interpolate_array_using_CIC_2d(nx,ny,xbdry,ybdry,tempin,newx,newy,sc_cnt,outvals);
               else if (interp == tsc)
                  interpolate_array_using_TSC_2d(nx,ny,xbdry,ybdry,tempin,newx,newy,sc_cnt,outvals);
               else
                  interpolate_array_using_M4p_2d(nx,ny,xbdry,ybdry,mask,tempin,newx,newy,sc_cnt,outvals);
               for (k=0; k<sc_cnt; k++) tempout[k][i][j] = outvals[k];
            } else {
               for (k=0; k<sc_cnt; k++) tempout[k][i][j] = tempin[k][i][j];
            }
         }
      }
   } else {
      fprintf(stderr,"Method-of-characteristics advection order %d unsupported.\n",order);
      fprintf(stderr," Try 1, 2, or 4.\n");
      fprintf(stderr," Quitting.\n");
      exit(0);
   }

   // prevent advecting colors
   if (!moveColors) {
      for (i=0;i<nx;i++) {
         for (j=0;j<ny;j++) {
            out[RR][i][j] = in[RR][i][j];
         }
      }
      for (i=0;i<nx;i++) {
         for (j=0;j<ny;j++) {
            out[GG][i][j] = in[GG][i][j];
         }
      }
      for (i=0;i<nx;i++) {
         for (j=0;j<ny;j++) {
            out[BB][i][j] = in[BB][i][j];
         }
      }
   }

   return(0);
}


/*
 * Interpolate velocity FROM the grid using the CIC method, but in 2
 * dimensions only, and with a vector-valued field
 *
 * loc is the point at which the vector is to be evaluated
 * zeta makes the 2d input field of vectors
 * out is the vector value interpolated from the field
 */
float interpolate_array_using_CIC_2d(int nx,int ny,int xbdry,int ybdry,
   float ***zeta,float px,float py,int numout,float out[MAX_SCALARS]) {

   int i,j,k,ii,ji,ir,jr;
   int si[2];           // start index
   float dx;
   float wt[2][2];
   float xfactor,yfactor;
   float mf;

   for (k=0; k<numout; k++) out[k] = 0.0;

   /* compute the lowest corner index of the cells that the scheme needs */
   // these may be negative
   si[0] = (int)(10+(nx-1)*(px)) - 10;
   si[1] = (int)(10+(nx-1)*(py)) - 10;
   // fprintf(stdout,"location is %g %g, start index is %d %d\n",px,py,si[0],si[1]);

   // make the list of CIC weights
   dx = fabs( (px)*(nx-1) - si[0] );
   wt[0][0] = 1.0-dx;
   dx = dx-1.0;
   wt[0][1] = 1.0+dx;

   // dx = fabs( (py)*(ny-1) - si[1] );
   dx = fabs( (py)*(nx-1) - si[1] );
   wt[1][0] = 1.0-dx;
   dx = dx-1.0;
   wt[1][1] = 1.0+dx;

   // ii,ji are the counters, i,j are the cell indexes (possibly negative)
   //   and ir,jr are the actual cells we take values from
   for (ii=0; ii<2; ii++) {
   for (ji=0; ji<2; ji++) {

      i = ii+si[0];
      j = ji+si[1];

      /* in order to properly reflect vorticity off walls, use xyzfactors */
      xfactor = 1.0;
      yfactor = 1.0;
      // but, when writing images, use additive reflection, i.e. no need for factors
      // also true when interpolating within scalar-type fields

      /* near the min or max x bounds */
      if (xbdry == PERIODIC) {
         ir = (i+(nx-1))%(nx-1);
      } else {
         ir = i;
         if (i<0) {
            ir = -i;
            // xfactor *= -1.0;
         } else if (i>(nx-1)) {
            ir = 2*(nx-1) - i;
            // xfactor *= -1.0;
         }
      }

      /* near the min or max y bounds */
      if (ybdry == PERIODIC) {
         jr = (j+(ny-1))%(ny-1);
      } else {
         jr = j;
         if (j<0) {
            jr = -j;
            // yfactor *= -1.0;
         } else if (j>(ny-1)) {
            jr = 2*(ny-1) - j;
            // yfactor *= -1.0;
         }
      }

      /* find the m4-factor for this cell's contribution */
      mf = wt[0][ii]*wt[1][ji]*xfactor*yfactor;

      /* apply them to the grid node in question */
      for (k=0; k<numout; k++) out[k] += zeta[k][ir][jr]*mf;

   }}

   /* all's well that ends well */
   return(0);
}


/*
 * Interpolate velocity FROM the grid using the TSC method, but in 2
 * dimensions only, and with a vector-valued field
 *
 * loc is the point at which the vector is to be evaluated
 * zeta makes the 2d input field of vectors
 * out is the vector value interpolated from the field
 */
float interpolate_array_using_TSC_2d(int nx,int ny,int xbdry,int ybdry,
   float ***zeta,float px,float py,int numout,float out[MAX_SCALARS]) {

   int i,j,k,ii,ji,ir,jr;
   int si[2];           // start index
   float dx;
   float m3[2][3];
   float xfactor,yfactor;
   float mf;

   for (k=0; k<numout; k++) out[k] = 0.0;

   /* compute the lowest corner index of the cells that the scheme needs */
   // these may be negative
   si[0] = (int)(30.5+(nx-1)*(px)) - 31;
   // si[1] = (int)(30+(ny-1)*(py)) - 31;
   si[1] = (int)(30.5+(nx-1)*(py)) - 31;
   //fprintf(stdout,"location is %g %g, start index is %d %d\n",px,py,si[0],si[1]);

   // make the list of M3 weights
   //dx = fabs( (px)*(nx-1) - si[0]) -2.;
   dx = si[0] - (px)*(nx-1);
   m3[0][0] = 0.5*pow(dx + 1.5,2);
   dx = dx+1.0;
   m3[0][1] = 0.75 - dx*dx;
   dx = dx+1.0;
   m3[0][2] = 0.5*pow(1.5 - dx,2);

   // dx = fabs( (py)*(ny-1) - si[1] );
   //dx = fabs( (py)*(nx-1) - si[1]) -2.;
   dx = si[1] - (py)*(nx-1);
   m3[1][0] = 0.5*pow(dx + 1.5,2);
   dx = dx+1.0;
   m3[1][1] = 0.75 - dx*dx;
   dx = dx+1.0;
   m3[1][2] = 0.5*pow(1.5 - dx,2);

   //fprintf(stdout,"  weights are %g %g %g and %g %g %g\n",m3[0][0],m3[0][1],m3[0][2],m3[1][0],m3[1][1],m3[1][2]);

   // ii,ji are the counters, i,j are the cell indexes (possibly negative)
   //   and ir,jr are the actual cells we take values from
   for (ii=0; ii<3; ii++) {
   for (ji=0; ji<3; ji++) {

      i = ii+si[0];
      j = ji+si[1];

      /* in order to properly reflect vorticity off walls, use xyzfactors */
      xfactor = 1.0;
      yfactor = 1.0;
      // but, when writing images, use additive reflection, i.e. no need for factors
      // also true when interpolating within scalar-type fields

      /* near the min or max x bounds */
      if (xbdry == PERIODIC) {
         ir = (i+(nx-1))%(nx-1);
      } else {
         ir = i;
         if (i<0) {
            ir = -i;
            // xfactor *= -1.0;
         } else if (i>(nx-1)) {
            ir = 2*(nx-1) - i;
            // xfactor *= -1.0;
         }
      }

      /* near the min or max y bounds */
      if (ybdry == PERIODIC) {
         jr = (j+(ny-1))%(ny-1);
      } else {
         jr = j;
         if (j<0) {
            jr = -j;
            // yfactor *= -1.0;
         } else if (j>(ny-1)) {
            jr = 2*(ny-1) - j;
            // yfactor *= -1.0;
         }
      }

      /* find the m3-factor for this cell's contribution */
      mf = m3[0][ii]*m3[1][ji]*xfactor*yfactor;

      /* apply them to the grid node in question */
      for (k=0; k<numout; k++) out[k] += zeta[k][ir][jr]*mf;

   }}

   /* all's well that ends well */
   return(0);
}


/*
 * Interpolate velocity FROM the grid using the M4' method, but in 2
 * dimensions only, and with a vector-valued field
 *
 * loc is the point at which the vector is to be evaluated
 * zeta makes the 2d input field of vectors
 * out is the vector value interpolated from the field
 */
float interpolate_array_using_M4p_2d(int nx,int ny,int xbdry,int ybdry,
      float **mask,float ***zeta,
      float px,float py,int numout,float out[MAX_SCALARS]) {

   int i,j,k,ii,ji,ir,jr,nm1;
   int si[2];           // start index
   float dx;
   float m4[2][4];
   float xfactor,yfactor;
   float mf;
   float mfsum = 0.;
   //const int debug = (px < 0.0 && py > 0.1 && py < 0.102);
   const int debug = FALSE;

   for (k=0; k<numout; k++) out[k] = 0.0;

   /* compute the lowest corner index of the cells that the scheme needs */
   // these may be negative
   if (nx > ny) {
      nm1 = nx-1;
   } else {
      nm1 = ny-1;
   }
   si[0] = (int)(20+(nm1)*(px)) - 21;
   si[1] = (int)(20+(nm1)*(py)) - 21;
   if (debug) fprintf(stdout,"location is %g %g, start index is %d %d\n",px,py,si[0],si[1]);

   // make the list of M4 weights
   dx = fabs( (px)*(nm1) - si[0] );
   m4[0][0] = 0.5*pow(2-dx,2)*(1-dx);
   dx = dx-1.0;
   //m4[0][1] = 1.0 - 2.5*dx*dx + 1.5*dx*dx*dx;
   m4[0][1] = 1.0 - dx*dx*(2.5 - 1.5*dx);
   dx = fabs(dx-1.0);
   //m4[0][2] = 1.0 - 2.5*dx*dx + 1.5*dx*dx*dx;
   m4[0][2] = 1.0 - dx*dx*(2.5 - 1.5*dx);
   dx = dx+1.0;
   m4[0][3] = 0.5*pow(2-dx,2)*(1-dx);

   // dx = fabs( (py)*(ny-1) - si[1] );
   dx = fabs( (py)*(nm1) - si[1] );
   m4[1][0] = 0.5*pow(2-dx,2)*(1-dx);
   dx = dx-1.0;
   //m4[1][1] = 1.0 - 2.5*dx*dx + 1.5*dx*dx*dx;
   m4[1][1] = 1.0 - dx*dx*(2.5 - 1.5*dx);
   dx = fabs(dx-1.0);
   //m4[1][2] = 1.0 - 2.5*dx*dx + 1.5*dx*dx*dx;
   m4[1][2] = 1.0 - dx*dx*(2.5 - 1.5*dx);
   dx = dx+1.0;
   m4[1][3] = 0.5*pow(2-dx,2)*(1-dx);

   if (debug) fprintf(stdout,"  weights are %g %g %g %g and %g %g %g %g\n",m4[0][0],m4[0][1],m4[0][2],m4[0][3],m4[1][0],m4[1][1],m4[1][2],m4[1][3]);

   // ii,ji are the counters, i,j are the cell indexes (possibly negative)
   //   and ir,jr are the actual cells we take values from
   for (ii=0; ii<4; ii++) {
   for (ji=0; ji<4; ji++) {

      i = ii+si[0];
      j = ji+si[1];

      /* in order to properly reflect vorticity off walls, use xyzfactors */
      //xfactor = 1.0;
      //yfactor = 1.0;
      // but, when writing images, use additive reflection, i.e. no need for factors
      // also true when interpolating within scalar-type fields

      /* near the min or max x bounds */
      ir = i;
      if (xbdry == PERIODIC) {
         ir = (i+(nx-1))%(nx-1);
      } else if (xbdry == OPEN) {
         if (i<0) {
            ir = 0;
         } else if (i>(nx-1)) {
            ir = nx-1;
         }
      } else {	// WALL
         if (i<0) {
            ir = -i;
            // xfactor *= -1.0;
         } else if (i>(nx-1)) {
            ir = 2*(nx-1) - i;
            // xfactor *= -1.0;
         }
      }

      /* near the min or max y bounds */
      jr = j;
      if (ybdry == PERIODIC) {
         jr = (j+(ny-1))%(ny-1);
      } else if (ybdry == OPEN) {
         if (j<0) {
            jr = 0;
         } else if (j>(ny-1)) {
            jr = ny-1;
         }
      } else {	// WALL
         if (j<0) {
            jr = -j;
            // yfactor *= -1.0;
         } else if (j>(ny-1)) {
            jr = 2*(ny-1) - j;
            // yfactor *= -1.0;
         }
      }

      /* find the m4-factor for this cell's contribution */
      mf = m4[0][ii] * m4[1][ji];
      if (mask != NULL) mf *= mask[ir][jr];
      mfsum += mf;

      if (debug) fprintf(stdout,"    ii,ji %d %d   ir,jr %d %d   weight %g\n",ii,ji,ir,jr,mf);

      /* apply them to the grid node in question */
      for (k=0; k<numout; k++) out[k] += zeta[k][ir][jr]*mf;

   }}

   // correct by the actual sum of the weights
   if (fabs(mfsum) > 0.) {
      for (k=0; k<numout; k++) out[k] /= mfsum;
   }
   if (debug) for (k=0; k<numout; k++) fprintf(stdout,"    k %d   val %g\n",k,out[k]);

   /* all's well that ends well */
   return(0);
}


/*
 * Interpolate velocity FROM the grid using the CIC method, but in 2
 * dimensions only, and with a vector-valued field
 *
 * loc is the point at which the vector is to be evaluated
 * zeta makes the 2d input field of vectors
 * out is the vector value interpolated from the field
 */
int interpolate_vel_using_CIC_2d(int nx,int ny,int xbdry,int ybdry,float **ua,float **va,float px,float py,float *u,float *v) {

   int i,j,ii,ji,ir,jr;
   int si[2];           // start index
   FLOAT dx;
   FLOAT wt[2][2];
   FLOAT xfactor,yfactor;
   FLOAT mf;
   FLOAT out = 0.0;

   /* compute the lowest corner index of the cells that the scheme needs */
   // these may be negative
   si[0] = (int)(10+(nx-1)*(px)) - 10;
   si[1] = (int)(10+(nx-1)*(py)) - 10;
   // fprintf(stdout,"location is %g %g, start index is %d %d\n",px,py,si[0],si[1]);

   // make the list of CIC weights
   dx = fabs( (px)*(nx-1) - si[0] );
   wt[0][0] = 1.0-dx;
   dx = dx-1.0;
   wt[0][1] = 1.0+dx;
                                                                                         
   dx = fabs( (py)*(nx-1) - si[1] );
   wt[1][0] = 1.0-dx;
   dx = dx-1.0;
   wt[1][1] = 1.0+dx;

   // fprintf(stdout,"  weights are %g %g %g %g and %g %g %g %g\n",m4[0][0],m4[0][1],m4[0][2],m4[0][3],m4[1][0],m4[1][1],m4[1][2],m4[1][3]);

   // ii,ji are the counters, i,j are the cell indexes (possibly negative)
   //   and ir,jr are the actual cells we take values from
   *(u) = 0.0;
   *(v) = 0.0;
   for (ii=0; ii<2; ii++) {
   for (ji=0; ji<2; ji++) {

      i = ii+si[0];
      j = ji+si[1];

      /* in order to properly reflect vorticity off walls, use xyzfactors */
      xfactor = 1.0;
      yfactor = 1.0;
      // but, when writing images, use additive reflection, i.e. no need for factors
      // also true when interpolating within scalar-type fields

      /* near the min or max x bounds */
      if (xbdry == PERIODIC) {
         ir = (i+(nx-1))%(nx-1);
      } else {
         ir = i;
         if (i==0) {
            // xfactor = 1.0;
         } else if (i<0) {
            ir = -i;
            // xfactor *= -1.0;
         } else if (i==(nx-1)) {
            // xfactor = 0.0;
         } else if (i>(nx-1)) {
            ir = 2*(nx-1) - i;
            // xfactor *= -1.0;
         }
      }

      /* near the min or max y bounds */
      if (ybdry == PERIODIC) {
         jr = (j+(ny-1))%(ny-1);
      } else {
         jr = j;
         if (j==0) {
            // yfactor = 1.0;
         } else if (j<0) {
            jr = -j;
            // yfactor *= -1.0;
         } else if (j==(ny-1)) {
            // yfactor = 0.0;
         } else if (j>(ny-1)) {
            jr = 2*(ny-1) - j;
            // yfactor *= -1.0;
         }
      }

      /* find the m4-factor for this cell's contribution */
      mf = xfactor*yfactor*wt[0][ii]*wt[1][ji];

      /* apply them to the grid node in question */
      *(u) += ua[ir][jr]*mf;
      *(v) += va[ir][jr]*mf;

   }}

   // fprintf(stdout,"      out is %g\n",out);

   /* all's well that ends well */
   return(0);
}


/*
 * Interpolate velocity FROM the grid using the M4' method, but in 2
 * dimensions only, and with a vector-valued field
 *
 * loc is the point at which the vector is to be evaluated
 * zeta makes the 2d input field of vectors
 * out is the vector value interpolated from the field
 */
int interpolate_vel_using_TSC_2d(int nx,int ny,int xbdry,int ybdry,float **ua,float **va,float px,float py,float *u,float *v) {

   int i,j,ii,ji,ir,jr;
   int si[2];           // start index
   FLOAT dx;
   FLOAT m3[2][3];
   FLOAT xfactor,yfactor;
   FLOAT mf;
   FLOAT out = 0.0;

   /* compute the lowest corner index of the cells that the scheme needs */
   // these may be negative
   si[0] = (int)(30.5+(nx-1)*(px)) - 31;
   // si[1] = (int)(30+(ny-1)*(py)) - 31;
   si[1] = (int)(30.5+(nx-1)*(py)) - 31;
   // fprintf(stdout,"location is %g %g, start index is %d %d\n",px,py,si[0],si[1]);

   // make the list of M4 weights
   dx = si[0] - (px)*(nx-1);
   m3[0][0] = 0.5*pow(dx + 1.5,2);
   dx = dx+1.0;
   m3[0][1] = 0.75 - dx*dx;
   dx = dx+1.0;
   m3[0][2] = 0.5*pow(1.5 - dx,2);

   dx = si[1] - (py)*(nx-1);
   m3[1][0] = 0.5*pow(dx + 1.5,2);
   dx = dx+1.0;
   m3[1][1] = 0.75 - dx*dx;
   dx = dx+1.0;
   m3[1][2] = 0.5*pow(1.5 - dx,2);

   // fprintf(stdout,"  weights are %g %g %g %g and %g %g %g %g\n",m4[0][0],m4[0][1],m4[0][2],m4[0][3],m4[1][0],m4[1][1],m4[1][2],m4[1][3]);

   // ii,ji are the counters, i,j are the cell indexes (possibly negative)
   //   and ir,jr are the actual cells we take values from
   *(u) = 0.0;
   *(v) = 0.0;
   for (ii=0; ii<3; ii++) {
   for (ji=0; ji<3; ji++) {

      i = ii+si[0];
      j = ji+si[1];

      /* in order to properly reflect vorticity off walls, use xyzfactors */
      xfactor = 1.0;
      yfactor = 1.0;
      // but, when writing images, use additive reflection, i.e. no need for factors
      // also true when interpolating within scalar-type fields

      /* near the min or max x bounds */
      if (xbdry == PERIODIC) {
         ir = (i+(nx-1))%(nx-1);
      } else {
         ir = i;
         if (i==0) {
            // xfactor = 1.0;
         } else if (i<0) {
            ir = -i;
            // xfactor *= -1.0;
         } else if (i==(nx-1)) {
            // xfactor = 0.0;
         } else if (i>(nx-1)) {
            ir = 2*(nx-1) - i;
            // xfactor *= -1.0;
         }
      }

      /* near the min or max y bounds */
      if (ybdry == PERIODIC) {
         jr = (j+(ny-1))%(ny-1);
      } else {
         jr = j;
         if (j==0) {
            // yfactor = 1.0;
         } else if (j<0) {
            jr = -j;
            // yfactor *= -1.0;
         } else if (j==(ny-1)) {
            // yfactor = 0.0;
         } else if (j>(ny-1)) {
            jr = 2*(ny-1) - j;
            // yfactor *= -1.0;
         }
      }

      /* find the m4-factor for this cell's contribution */
      mf = xfactor*yfactor*m3[0][ii]*m3[1][ji];

      /* apply them to the grid node in question */
      *(u) += ua[ir][jr]*mf;
      *(v) += va[ir][jr]*mf;

   }}

   // fprintf(stdout,"      out is %g\n",out);

   /* all's well that ends well */
   return(0);
}


/*
 * Interpolate velocity FROM the grid using the M4' method, but in 2
 * dimensions only, and with a vector-valued field
 *
 * loc is the point at which the vector is to be evaluated
 * zeta makes the 2d input field of vectors
 * out is the vector value interpolated from the field
 */
int interpolate_vel_using_M4p_2d(int nx,int ny,int xbdry,int ybdry,
      float **mask,float **ua,float **va,
      float px,float py,float *u,float *v) {

   int i,j,ii,ji,ir,jr,nm1;
   int si[2];           // start index
   FLOAT dx;
   FLOAT m4[2][4];
   FLOAT xfactor,yfactor;
   FLOAT mf;
   FLOAT out = 0.0;
   FLOAT mfsum = 0.;

   /* compute the lowest corner index of the cells that the scheme needs */
   // these may be negative
   if (nx > ny) {
      nm1 = nx-1;
   } else {
      nm1 = ny-1;
   }
   si[0] = (int)(20+(nm1)*(px)) - 21;
   si[1] = (int)(20+(nm1)*(py)) - 21;
   // fprintf(stdout,"location is %g %g, start index is %d %d\n",px,py,si[0],si[1]);

   // make the list of M4 weights
   // first in the x direction
   dx = fabs( (px)*(nm1) - si[0] );
   m4[0][0] = 0.5*pow(2.0-dx,2)*(1.0-dx);
   dx = dx-1.0;
   //m4[0][1] = 1.0 - 2.5*dx*dx + 1.5*dx*dx*dx;
   m4[0][1] = 1.0 + dx*dx*(1.5*dx - 2.5);
   dx = fabs(dx-1.0);
   //m4[0][2] = 1.0 - 2.5*dx*dx + 1.5*dx*dx*dx;
   m4[0][2] = 1.0 + dx*dx*(1.5*dx - 2.5);
   dx = dx+1.0;
   m4[0][3] = 0.5*pow(2.0-dx,2)*(1.0-dx);

   // then in the y direction
   dx = fabs( (py)*(nm1) - si[1] );
   m4[1][0] = 0.5*pow(2.0-dx,2)*(1.0-dx);
   dx = dx-1.0;
   //m4[1][1] = 1.0 - 2.5*dx*dx + 1.5*dx*dx*dx;
   m4[1][1] = 1.0 + dx*dx*(1.5*dx - 2.5);
   dx = fabs(dx-1.0);
   //m4[1][2] = 1.0 - 2.5*dx*dx + 1.5*dx*dx*dx;
   m4[1][2] = 1.0 + dx*dx*(1.5*dx - 2.5);
   dx = dx+1.0;
   m4[1][3] = 0.5*pow(2.0-dx,2)*(1.0-dx);

   // fprintf(stdout,"  weights are %g %g %g %g and %g %g %g %g\n",m4[0][0],m4[0][1],m4[0][2],m4[0][3],m4[1][0],m4[1][1],m4[1][2],m4[1][3]);

   // ii,ji are the counters, i,j are the cell indexes (possibly negative)
   //   and ir,jr are the actual cells we take values from
   *(u) = 0.0;
   *(v) = 0.0;
   for (ii=0; ii<4; ii++) {
   for (ji=0; ji<4; ji++) {

      i = ii+si[0];
      j = ji+si[1];

      /* in order to properly reflect vorticity off walls, use xyzfactors */
      //xfactor = 1.0;
      //yfactor = 1.0;
      // but, when writing images, use additive reflection, i.e. no need for factors
      // also true when interpolating within scalar-type fields

      /* near the min or max x bounds */
      ir = i;
      if (xbdry == PERIODIC) {
         ir = (i+(nx-1))%(nx-1);
      } else if (xbdry == OPEN) {
         if (i<0) {
            ir = 0;
         } else if (i>(nx-1)) {
            ir = nx-1;
         }
      } else {
         if (i==0) {
            // xfactor = 1.0;
         } else if (i<0) {
            ir = -i;
            // xfactor *= -1.0;
         } else if (i==(nx-1)) {
            // xfactor = 0.0;
         } else if (i>(nx-1)) {
            ir = 2*(nx-1) - i;
            // xfactor *= -1.0;
         }
      }

      /* near the min or max y bounds */
      jr = j;
      if (ybdry == PERIODIC) {
         jr = (j+(ny-1))%(ny-1);
      } else if (ybdry == OPEN) {
         if (j<0) {
            jr = 0;
         } else if (j>(ny-1)) {
            jr = ny-1;
         }
      } else {
         if (j==0) {
            // yfactor = 1.0;
         } else if (j<0) {
            jr = -j;
            // yfactor *= -1.0;
         } else if (j==(ny-1)) {
            // yfactor = 0.0;
         } else if (j>(ny-1)) {
            jr = 2*(ny-1) - j;
            // yfactor *= -1.0;
         }
      }

      /* find the m4-factor for this cell's contribution */
      mf = m4[0][ii] * m4[1][ji];
      if (mask != NULL) mf *= mask[ir][jr];
      mfsum += mf;

      /* apply them to the grid node in question */
      *(u) += ua[ir][jr]*mf;
      *(v) += va[ir][jr]*mf;

   }}

   // correct by the actual sum of the weights
   if (fabs(mfsum) > 0.) {
      *(u) /= mfsum;
      *(v) /= mfsum;
   }

   // fprintf(stdout,"      out is %g\n",out);

   /* all's well that ends well */
   return(0);
}


