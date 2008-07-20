/*
 * VIC-MOC - libvicmoc.c - the library for the 2D and 3D routines
 *
 * Copyright 2004-8 Mark J. Stock mstock@umich.edu
 */

#include "vicmoc.h"
#ifdef _OPENMP
#include <omp.h>
#endif


int add_sharp_circular_blob(int,int,int,int,float**,float,float,float,float);
int add_smooth_circular_blob(int,int,int,int,float**,float,float,float,float);
int add_smooth_circular_blob_3d(int,int,int,int,int,int,int,float***,float,float,float,float);
int add_smooth_spherical_blob(int,int,int,int,int,int,float***,float,float,float,float,float);
int add_singular_blob_3d(int,int,int,int,int,int,float***,float,float,float,float,float);
int add_cube_3d(int,int,int,int,int,int,float***,float,float,float,float,float);
int add_vortex_ring_3d(int,int,int,int,int,int,float***,float***,float***,float,float,float,float,float,float,float,float,float);
int explicit_particle_move_3d(int,int,int,int,int,int,float****,float,float,int,float**,float**);
float step_forward_2d(int,int,int,int,int,int,float*,float***,float***,float***,int,float**,int,int*,float*,int,float*,float,float);
float step_forward_3d(int,int,int,int,int,int,int,int,float****,float****,int,int*,float*,float,float);

int create_baroclinic_vorticity_2d(int,int,int,int,float**,float**,float,float*,float);
int create_baroclinic_vorticity_3d(int,int,int,int,int,int,float***,float***,float***,float***,float,float);
int create_vorticity_stretch_3d(int,int,int,int,int,int,float***,float***,float***,float***,float***,float***,float);
int create_boundary_vorticity_2d(int,int,int,int,float**,float**,float**,float);
int create_boundary_vorticity_3d(int,int,int,int,int,int,float***,float***,float***,float***,float***,float***,float);
float diffuse_scalar_2d(int,int,int,int,float**,float**,float,int,float**,float);
float diffuse_scalar_3d(int,int,int,int,int,int,float***,float***,float,float);
int find_vels_2d(int,int,int,int,int,int,float*,float**,float**,float**,int,float**);
int find_vels_3d(int,int,int,int,int,int,int,float***,float***,float***,float***,float***,float***);
int make_solenoidal_3d(int,int,int,int,int,int,float***,float***,float***);
int find_gradient_of_scalar_2d(int,int,int,int,float**,float**,float,float**,float);
int find_gradient_of_scalar_2nd_2d(int,int,int,int,float**,float**,float,float**,float);
int find_gradient_of_scalar_3d(int,int,int,int,int,int,float***,float****);
int find_curl_of_vel_2d(int,int,int,int,float**,float**,float**);
int moc_advect_2d(int,int,int,int,float**,float**,float***,float***,int,float);
int moc_advect_3d(int,int,int,int,int,int,float****,float****,float****,int,float);
float compute_divergence_3d(int,int,int,int,int,int,float***,float***,float***,float***);
int compute_curl_3d(int,int,int,int,int,int,float***,float***,float***,float***,float***,float***);
float find_energy_3d(int,int,int,int,int,int,float***,float***,float***);
float find_biot_savart_u (float,float,int,int,float,float,float**);
float find_biot_savart_v (float,float,int,int,float,float,float**);
int find_biot_savart (float,float,int,int,float,float,float**,float*);
int find_open_boundary_psi (int,int,float**,float,float*,float**);
float interpolate_using_area_2d(int,int,int,int,float**,float,float);
float interpolate_array_using_CIC_2d(int,int,int,int,float***,float,float,int,float*);
int interpolate_vel_using_CIC_2d(int,int,int,int,float**,float**,float,float,float*,float*);
float interpolate_array_using_TSC_2d(int,int,int,int,float***,float,float,int,float*);
int interpolate_vel_using_TSC_2d(int,int,int,int,float**,float**,float,float,float*,float*);
float interpolate_using_M4p_2d(int,int,int,int,float**,float,float);
float interpolate_array_using_M4p_2d(int,int,int,int,float***,float,float,int,float*);
int interpolate_vel_using_M4p_2d(int,int,int,int,float**,float**,float,float,float*,float*);
float interpolate_array_using_CIC_3d(int,int,int,int,int,int,float****,float,float,float,int,float*);
float interpolate_array_using_M4p_3d(int,int,int,int,int,int,float****,float,float,float,int,float*);

float* allocate_1d_array_f(int);
int free_1d_array_f(float*);
float** allocate_2d_array_f(int,int);
int free_2d_array_f(float**);
float*** allocate_3d_array_f(int,int,int);
int free_3d_array_f(float***);
double gaussian(double,double);

int copy_2d_field(float**,float**,int,int);


/*
 * For whatever reason, add a chunk of scalar stuff
 */
int add_sharp_circular_blob(int nx,int ny,int xbdry,int ybdry,float **addto,float xpos,float ypos,float rad,float val) {

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
 * add the same thing, but with wy vorticity
 */
int add_smooth_circular_blob_3d(int nx,int ny,int nz,int xbdry,int ybdry,int zbdry,int axis,float ***addto,float xpos,float ypos,float rad,float val) {

   int ix,iy,iz;
   float px,py,pz;
   float rad_inner = 0.8*rad;
   float dist;
   float temp;

   if (axis==0) {
     // local x is y, local y is z
     for (iy=0; iy<ny; iy++) {
       py = (float)iy/(float)(ny-1);
       for (iz=0; iz<nz; iz++) {
         pz = (float)iz/(float)(nz-1);
         dist = sqrt(pow(py-xpos,2)+pow(pz-ypos,2));
         if (dist < rad_inner) {
           for (ix=0; ix<nx; ix++) {
             addto[ix][iy][iz] += val;
           }
         } else if (dist < rad) {
           temp = val*(0.5+0.5*cos(M_PI*(dist-rad_inner)/(rad-rad_inner)));
           for (ix=0; ix<nx; ix++) {
             addto[ix][iy][iz] += temp;
           }
         }
       }
     }
   } else if (axis==1) {
     // local x is x, local y is z
     for (ix=0; ix<nx; ix++) {
       px = (float)ix/(float)(nx-1);
       for (iz=0; iz<nz; iz++) {
         py = (float)iz/(float)(nz-1);
         dist = sqrt(pow(px-xpos,2)+pow(py-ypos,2));
         if (dist < rad_inner) {
           for (iy=0; iy<ny; iy++) {
             addto[ix][iy][iz] += val;
           }
         } else if (dist < rad) {
           temp = val*(0.5+0.5*cos(M_PI*(dist-rad_inner)/(rad-rad_inner)));
           for (iy=0; iy<ny; iy++) {
             addto[ix][iy][iz] += temp;
           }
         }
       }
     }
   } else if (axis==2) {
     // local x is x, local y is y
     for (ix=0; ix<nx; ix++) {
       px = (float)ix/(float)(nx-1);
       for (iy=0; iy<ny; iy++) {
         py = (float)iy/(float)(ny-1);
         dist = sqrt(pow(px-xpos,2)+pow(py-ypos,2));
         if (dist < rad_inner) {
           for (iz=0; iz<nz; iz++) {
             addto[ix][iy][iz] += val;
           }
         } else if (dist < rad) {
           temp = val*(0.5+0.5*cos(M_PI*(dist-rad_inner)/(rad-rad_inner)));
           for (iz=0; iz<nz; iz++) {
             addto[ix][iy][iz] += temp;
           }
         }
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
 * Add a sharp (high peak) Gaussian scalar distribution, hopefully like the
 * temp field after a nuclear explosion (or any explosion, really)
 *
 * use Rosenhead/Krasny desingularization
 */
int add_singular_blob_3d(int nx,int ny,int nz,int xbdry,int ybdry,int zbdry,
   float ***addto,float xpos,float ypos,float zpos,float rad,float val) {

   int ix,iy,iz;
   float px,py,pz;
   float epsilon = 1.e-5;
   float offset = 1./(fabs(pow(rad,3))+epsilon);
   float dist;

   for (ix=0; ix<nx; ix++) {
      px = (float)ix/(float)(nx-1);
      for (iy=0; iy<ny; iy++) {
         py = (float)iy/(float)(nx-1);
         for (iz=0; iz<nz; iz++) {
            pz = (float)iz/(float)(nx-1);
            dist = sqrt(pow(px-xpos,2)+pow(py-ypos,2)+pow(pz-zpos,2));
            if (dist < rad) {
               addto[ix][iy][iz] += val*(1./(fabs(pow(dist,3))+epsilon)-offset);
            }
         }
      }
   }

   return(0);
}


/*
 * Add a cube of scalar
 */
int add_cube_3d(int nx,int ny,int nz,int xbdry,int ybdry,int zbdry,
   float ***addto,float xpos,float ypos,float zpos,float size,float val) {

   int ix,iy,iz;
   int sx,sy,sz;
   int fx,fy,fz;

   sx = (int)((float)(xpos)*(float)(nx-1));
   sy = (int)((float)(ypos)*(float)(nx-1));
   sz = (int)((float)(zpos)*(float)(nx-1));
   fx = (int)((float)(xpos+size)*(float)(nx-1));
   fy = (int)((float)(ypos+size)*(float)(nx-1));
   fz = (int)((float)(zpos+size)*(float)(nx-1));
   if (sx<0) sx = 0;
   if (sy<0) sy = 0;
   if (sz<0) sz = 0;
   if (fx>nx-1) fx = nx-1;
   if (fy>ny-1) fy = ny-1;
   if (fz>nz-1) fz = nz-1;

   for (ix=sx; ix<fx+1; ix++) {
      for (iy=sy; iy<fy+1; iy++) {
         for (iz=sz; iz<fz+1; iz++) {
            addto[ix][iy][iz] += val;
         }
      }
   }

   return(0);
}


/*
 * add the same thing, but with wy vorticity
 */
int add_vortex_ring_3d(int nx,int ny,int nz,int xbdry,int ybdry,int zbdry,
    float ***wx,float ***wy,float ***wz,float xloc,float yloc,float zloc,
    float xn,float yn,float zn,float majrad,float minrad,float str) {

   int ix,iy,iz;
   float px,py,pz;
   float rx,ry,rz,lx,ly,lz,zl,rl;
   float dist,temp,len,cx,cy,cz;

   // just to be safe, normalize the normal
   len = sqrt(xn*xn+yn*yn+zn*zn);
   xn /= len;
   yn /= len;
   zn /= len;

   // just loop over all positions
   for (ix=0; ix<nx; ix++) {
     px = (float)ix/(float)(nx-1);
     for (iy=0; iy<ny; iy++) {
       py = (float)iy/(float)(ny-1);
       for (iz=0; iz<nz; iz++) {
         pz = (float)iz/(float)(nz-1);

         // now, figure out where we are!
         //printf("point at %g %g %g\n",px,py,pz);

         // relative vector to the ring center
         rx = px-xloc;
         ry = py-yloc;
         rz = pz-zloc;
         //printf("  rel pos %g %g %g\n",rx,ry,rz);

         // distance from the plane
         zl = rx*xn + ry*yn + rz*zn;
         //printf("  dist from plane %g\n",zl);

         // location on the plane
         lx = rx-zl*xn;
         ly = ry-zl*yn;
         lz = rz-zl*zn;
         //printf("  inplane pos %g %g %g\n",lx,ly,lz);

         // distance (in the plane) from the ring axis
         rl = sqrt(lx*lx+ly*ly+lz*lz);
         //printf("  dist from axis %g\n",rl);

         // distance from the ring core
         dist = sqrt(zl*zl+pow(rl-majrad,2));
         //printf("  dist from core %g\n",dist);

         if (dist < minrad) {
           // strength of this point
           temp = str*(0.5+0.5*cos(M_PI*dist/minrad));

           // normal cross local (in-plane) vector is direction
           cx = yn*lz - zn*ly;
           cy = zn*lx - xn*lz;
           cz = xn*ly - yn*lx;
           //printf("  str dir %g %g %g\n",cx,cy,cz);
           //exit(0);

           // normalize this one, too
           dist = sqrt(cx*cx + cy*cy + cz*cz);
           cx /= dist;
           cy /= dist;
           cz /= dist;

           // apply it to the vorticities
           wx[ix][iy][iz] += cx*temp;
           wy[ix][iy][iz] += cy*temp;
           wz[ix][iy][iz] += cz*temp;
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
float step_forward_2d (int step, int isStam,
    int nx, int ny, int xbdry, int ybdry, float *freestream,
    float ***u, float ***a, float ***t,
    int use_MASK,float **mask,
    int sc_cnt, int sc_type[MAX_SCALARS], float sc_diffus[MAX_SCALARS],
    int gravtype, float *grav,
    float dt, float bn) {

  int i,j,k;
  int sfi = -1;
  int ufi = -1;
  int vfi = -1;
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
    for (i=0; i<MAX_SCALARS; i++) if (sc_type[i] == XV) ufi = i;
    (void) copy_2d_field(u[XV], a[ufi], nx, ny);
    for (i=0; i<MAX_SCALARS; i++) if (sc_type[i] == YV) vfi = i;
    (void) copy_2d_field(u[YV], a[vfi], nx, ny);

    //sprintf(outfileroot,"u_%04d",step+10);
    //write_output(outfileroot,nx,ny,a[ufi],-0.25,0.5,1);
    //sprintf(outfileroot,"v_%04d",step+10);
    //write_output(outfileroot,nx,ny,a[vfi],-0.25,0.5,1);

  } else {

    // This uses the new, unpublished method by Mark Stock whereupon the
    //  new *vorticity* field is updated in the same manner as Stam's
    //  velocity field.

    // create new vorticity based on density differences (gradient of scalar)
    for (i=1; i<MAX_SCALARS; i++) if (sc_type[i] == SF) sfi = i;
    if (sfi > -1)
      create_baroclinic_vorticity_2d (nx,ny,xbdry,ybdry,a[sfi],a[W2],dt,grav,bn);

    // create new vorticity due to viscous boundaries, etc,
    // sets 'a' vorticity at boundaries
    create_boundary_vorticity_2d (nx,ny,xbdry,ybdry,u[XV],u[YV],a[W2],dt);

  }

  // diffuse each scalar individually (a[] becomes t[])
  for (k=0; k<sc_cnt; k++) {
    if (TRUE) {
      coeff = diffuse_scalar_2d (nx,ny,xbdry,ybdry,a[k],t[k],sc_diffus[k],use_MASK,mask,dt);
    } else {
      // if we don't diffuse it, we still need to copy the values
      for (i=0;i<nx;i++) for (j=0;j<ny;j++) {
        t[k][i][j] = a[k][i][j];
      }
      coeff = 1.0;
    }
  }

  if (isStam) {
    //sprintf(outfileroot,"u_%04d",step+20);
    //write_output(outfileroot,nx,ny,t[ufi],-0.25,0.5,1);
    //sprintf(outfileroot,"v_%04d",step+20);
    //write_output(outfileroot,nx,ny,t[vfi],-0.25,0.5,1);
  }

  // again, split on method
  if (isStam) {

    // must copy 't' back to 'u' so that it can be used in moc advection
    (void) copy_2d_field(t[ufi], u[XV], nx, ny);
    (void) copy_2d_field(t[vfi], u[YV], nx, ny);

    //sprintf(outfileroot,"u_%04d",step+30);
    //write_output(outfileroot,nx,ny,u[XV],-0.25,0.5,1);
    //sprintf(outfileroot,"v_%04d",step+30);
    //write_output(outfileroot,nx,ny,u[YV],-0.25,0.5,1);

  }

  // project forward to find the new fields (both)
  fprintf(stderr,"  now in moc_advect_2d\n"); fflush(stderr);
  moc_advect_2d (nx,ny,xbdry,ybdry,u[XV],u[YV],t,a,sc_cnt,dt);

  if (isStam) {
    //sprintf(outfileroot,"u_%04d",step+40);
    //write_output(outfileroot,nx,ny,a[ufi],-0.25,0.5,1);
    //sprintf(outfileroot,"v_%04d",step+40);
    //write_output(outfileroot,nx,ny,a[vfi],-0.25,0.5,1);
  }

  // split on method for advection step
  if (isStam) {

    // apply no-slip wall boundary conditions
    if (ybdry == WALL) {
      j = ny-1;
      for (i=0;i<nx;i++) {
         a[ufi][i][j] = 0.;
         a[vfi][i][j] = 0.;
      }
      j = 0;
      for (i=0;i<nx;i++) {
         a[ufi][i][j] = 0.;
         a[vfi][i][j] = 0.;
      }
    }
    if (xbdry == WALL) {
      i = 0;
      for (j=0;j<ny;j++) {
         a[ufi][i][j] = 0.;
         a[vfi][i][j] = 0.;
      }
      i = nx-1;
      for (j=0;j<ny;j++) {
         a[ufi][i][j] = 0.;
         a[vfi][i][j] = 0.;
      }
    }

    // already used 'u' to find 'a' field from 't' field

    // calculate the vorticity field from 'a' field
    find_curl_of_vel_2d (nx,ny,xbdry,ybdry,a[ufi],a[vfi],a[W2]);

    // project forward to find the new velocity field
    // project 'a' velocity field onto a divergence-free vel field 'u'

    // make sure that the true velocities are in 'u' before going back

  } else {

    // create new vorticity due to viscous boundaries, etc,
    // sets 'a' vorticity at boundaries
    create_boundary_vorticity_2d (nx,ny,xbdry,ybdry,u[XV],u[YV],a[W2],dt);

  }

  // find the velocity field from the vorticity field
  fprintf(stderr,"  now in find_vels_2d\n"); fflush(stderr);
  find_vels_2d (step,isStam,nx,ny,xbdry,ybdry,freestream,u[XV],u[YV],a[W2],use_MASK,mask);

  // calculate effective Reynolds number
  Re = pow((float)(nx-1),2)*dt/coeff;

  // return Re, why? I dunno.
  return(Re);
}


/*
 * Procedure to take one step forward, given only vorticity
 * and scalar fraction - for 3D
 */
float step_forward_3d (int step,int isStam,
   int nx,int ny,int nz,int xbdry,int ybdry,int zbdry,
   float ****u,float ****a,int sc_cnt,
   int sc_type[MAX_SCALARS],float sc_diffus[MAX_SCALARS],
   float dt,float bn) {

   int i,j,k,l;
   int sfi = -1;
   float Re,coeff;
   static float ***t[MAX_SCALARS];
   static int set_temp = FALSE;

   // die if isStam
   if (isStam) {
     printf("ERROR (step_forward_3d): Stam forward advection not programmed\n");
     exit(0);
   }

   // malloc temporary arrays
   if (!set_temp) {
      for (i=0; i<sc_cnt; i++) {
         t[i] = allocate_3d_array_f(nx,ny,nz);
      }
      set_temp = TRUE;
   }

   // ----------------------------
   // find the velocity field from the vorticity field
   find_vels_3d(step,nx,ny,nz,xbdry,ybdry,zbdry,u[XV],u[YV],u[ZV],a[WX],a[WY],a[WZ]);
   //printf(" after find_vels %g %g %g\n",a[WY][7][0][8],a[WY][7][8][8],a[WY][7][16][8]);
   //printf("            vels %g %g %g\n",u[XV][7][0][8],u[YV][7][8][8],u[ZV][7][16][8]);
   //printf(" after find_vels %g %g %g\n",a[WY][8][0][8],a[WY][8][8][8],a[WY][8][16][8]);
   //printf("after find_vels %g %g %g\n",u[XV][8][16][8],u[YV][8][16][8],u[ZV][8][16][8]);
   //printf("           vort %g %g %g\n",a[WX][8][16][8],a[WY][8][16][8],a[WZ][8][16][8]);

   // project forward to find the new vorticities
   moc_advect_3d(nx,ny,nz,xbdry,ybdry,zbdry,u,a,t,sc_cnt,dt);
   //printf(" after moc advect %g %g %g\n",t[WY][8][8][0],t[WY][8][8][1],t[WY][8][8][2]);
   //printf(" after moc_advect_3d %g %g %g\n",t[WY][7][0][8],t[WY][7][8][8],t[WY][7][16][8]);
   //printf(" after moc_advect_3d %g %g %g\n",t[WY][8][0][8],t[WY][8][8][8],t[WY][8][16][8]);

   // moc does not account for stretch, do it here
   // NO! moc automatically accounts for stretch! think of 
   // vortex stretch as similar to evolution equation for
   // differential line element---both have stretching term.
   // moc routine automatically scales vorticity.
   // WRONG: MOC does not account for stretch! Must do this:
   // Yes, projecting onto solenoidal doesn't create stretch either!
   create_vorticity_stretch_3d(nx,ny,nz,xbdry,ybdry,zbdry,u[XV],u[YV],u[ZV],t[WX],t[WY],t[WZ],dt);
   // Maybe a better idea is solving a poisson equation to project
   // these divergent vorticity vectors onto a divergence-free field?
   make_solenoidal_3d(nx,ny,nz,xbdry,ybdry,zbdry,t[WX],t[WY],t[WZ]);
   //printf(" after stretch %g %g %g\n",t[WY][7][0][8],t[WY][7][8][8],t[WY][7][16][8]);
   //printf(" after stretch %g %g %g\n",t[WY][8][0][8],t[WY][8][8][8],t[WY][8][16][8]);

   // create new vorticity due to viscous boundaries, etc, modifies temporary vorticity
   create_boundary_vorticity_3d(nx,ny,nz,xbdry,ybdry,zbdry,u[XV],u[YV],u[ZV],t[WX],t[WY],t[WZ],dt);
   //printf(" after boundary vort %g %g %g\n",t[WY][8][8][0],t[WY][8][8][1],t[WY][8][8][2]);
   //printf(" after boundary vort %g %g %g\n",t[WY][8][0][8],t[WY][8][8][8],t[WY][8][16][8]);
   //printf(" after boundary vort %g %g %g\n",t[WY][7][0][8],t[WY][7][8][8],t[WY][7][16][8]);

   // create new vorticity based on density differences (gradient of scalar)
   for (i=1; i<MAX_SCALARS; i++) if (sc_type[i] == SF) sfi = i;
   if (sfi > -1)
      create_baroclinic_vorticity_3d(nx,ny,nz,xbdry,ybdry,zbdry,a[sfi],t[WX],t[WY],t[WZ],dt,bn);

   // diffuse each scalar individually (t[] becomes a[] again)
   for (l=0; l<sc_cnt; l++) {
      if (TRUE) {
         coeff = diffuse_scalar_3d(nx,ny,nz,xbdry,ybdry,zbdry,t[l],a[l],sc_diffus[l],dt);
      } else {
         // if we don't diffuse it, we still need to copy the values
         for (i=0;i<nx;i++) for (j=0;j<ny;j++) for (k=0;k<nz;k++) {
            a[l][i][j][k] = t[l][i][j][k];
         }
         coeff = 1.0;
      }
   }
   //printf(" after diffuse %g %g %g\n",t[WY][8][8][0],t[WY][8][8][1],t[WY][8][8][2]);
   //printf(" after diffuse %g %g %g\n",a[WY][8][0][8],a[WY][8][8][8],a[WY][8][16][8]);
   //printf(" after diffuse %g %g %g\n",a[WY][7][0][8],a[WY][7][8][8],a[WY][7][16][8]);

   // calculate effective Reynolds number
   Re = pow((float)(nx-1),2)*dt/coeff;

   // return 0, why? I dunno.
   return(Re);
}


/*
 * create vorticity at areas of scalar gradient
 */
int create_baroclinic_vorticity_2d (int nx,int ny,int xbdry,int ybdry,
      float **scalar,float **vort,float dt,float *g,float bn) {

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
   find_gradient_of_scalar_2nd_2d (nx,ny,xbdry,ybdry,scalar,grad[0],bn,grad[1],bn);

   // then, cross with the gravity vector
   for (i=0;i<nx;i++) {
      for (j=0;j<ny;j++) {
         vort[i][j] += dt*(grad[0][i][j]*g[1] - grad[1][i][j]*g[0]);
      }
   }

   // nah, just keep it. it's OK. it needs to diffuse, anyway
   return(0);

   // new way: treat boundaries properly!
   if (ybdry == WALL) {
      // correct top and bottom by taking the effect away
      j = ny-1;
      for (i=0;i<nx;i++) vort[i][j] -= dt*(grad[0][i][j]*g[1] - grad[1][i][j]*g[0]);
      j = 0;
      for (i=0;i<nx;i++) vort[i][j] -= dt*(grad[0][i][j]*g[1] - grad[1][i][j]*g[0]);
   }
   if (xbdry == WALL) {
      // correct left and right by taking the effect away
      i = 0;
      for (j=0;j<ny;j++) vort[i][j] -= dt*(grad[0][i][j]*g[1] - grad[1][i][j]*g[0]);
      i = nx-1;
      for (j=0;j<ny;j++) vort[i][j] -= dt*(grad[0][i][j]*g[1] - grad[1][i][j]*g[0]);
   }
   if (xbdry == WALL && ybdry == WALL) {
      // re-correct by adding it back on (we took it off twice!
      vort[0][0] += dt*(grad[0][0][0]*g[1] - grad[1][0][0]*g[0]);
      vort[0][ny-1] += dt*(grad[0][0][ny-1]*g[1] - grad[1][0][ny-1]*g[0]);
      vort[nx-1][0] += dt*(grad[0][nx-1][0]*g[1] - grad[1][nx-1][0]*g[0]);
      vort[nx-1][ny-1] += dt*(grad[0][nx-1][ny-1]*g[1] - grad[1][nx-1][ny-1]*g[0]);
   }

   return(0);
}


/*
 * create vorticity at areas of scalar gradient
 */
int create_baroclinic_vorticity_3d(int nx,int ny,int nz,
   int xbdry,int ybdry,int zbdry,
   float ***scalar,float ***wx,float ***wy,float ***wz,
   float dt, float bn) {

   int i,j,k;
   float g[3] = {0.0,0.0,-1.0};
   static float ***grad[3];
   static int set_grad = FALSE;

   // fprintf(stderr,"in create_baroclinic_vorticity_3d\n"); fflush(stderr);

   // if no boussinesq number (bn ~ 0), skip this!
   if (fabs(bn) < EPSILON) return(0);

   for (i=0; i<3; i++) g[i] = bn*g[i];

   // allocate memory for the gradient
   if (!set_grad) {
      grad[0] = allocate_3d_array_f(nx,ny,nz);
      grad[1] = allocate_3d_array_f(nx,ny,nz);
      grad[2] = allocate_3d_array_f(nx,ny,nz);
      set_grad = TRUE;
   }

   // first, find the scalar gradient
   find_gradient_of_scalar_3d (nx,ny,nz,xbdry,ybdry,zbdry,scalar,grad);

   // then, cross with the gravity vector (what about wall boundaries?)
   for (i=0;i<nx;i++) {
      for (j=0;j<ny;j++) {
         for (k=0;k<nz;k++) {
            wx[i][j][k] += dt*(grad[1][i][j][k]*g[2] - grad[2][i][j][k]*g[1]);
            wy[i][j][k] += dt*(grad[2][i][j][k]*g[0] - grad[0][i][j][k]*g[2]);
            wz[i][j][k] += dt*(grad[0][i][j][k]*g[1] - grad[1][i][j][k]*g[0]);
         }
      }
      // fprintf(stderr,"%g %g %g\n",wy[i][32][0],wy[i][32][1],wy[i][32][2]);
   }

   // nah, just keep it. it's OK. it needs to diffuse, anyway
   return(0);
}


/*
 * create vorticity due to vortex stretching
 *
 * first version: use velgrad at the old time but new locations
 * better: use velgrad at old locations, old vorticities?
 */
int create_vorticity_stretch_3d(int nx,int ny,int nz,
   int xbdry,int ybdry,int zbdry,
   float ***u,float ***v,float ***w,
   float ***wx,float ***wy,float ***wz,float dt) {

   int i,j,k;
   float wcgu[3];
   static float ***grad[3][3];
   static int set_grad = FALSE;

   // fprintf(stderr,"in create_vorticity_stretch_3d\n"); fflush(stderr);

   // allocate memory for the gradient
   if (!set_grad) {
      grad[0][0] = allocate_3d_array_f(nx,ny,nz);
      grad[0][1] = allocate_3d_array_f(nx,ny,nz);
      grad[0][2] = allocate_3d_array_f(nx,ny,nz);
      grad[1][0] = allocate_3d_array_f(nx,ny,nz);
      grad[1][1] = allocate_3d_array_f(nx,ny,nz);
      grad[1][2] = allocate_3d_array_f(nx,ny,nz);
      grad[2][0] = allocate_3d_array_f(nx,ny,nz);
      grad[2][1] = allocate_3d_array_f(nx,ny,nz);
      grad[2][2] = allocate_3d_array_f(nx,ny,nz);
      set_grad = TRUE;
   }

   // first, find the scalar gradient of the velocity
   find_gradient_of_scalar_3d (nx,ny,nz,xbdry,ybdry,zbdry,u,grad[0]);
   find_gradient_of_scalar_3d (nx,ny,nz,xbdry,ybdry,zbdry,v,grad[1]);
   find_gradient_of_scalar_3d (nx,ny,nz,xbdry,ybdry,zbdry,w,grad[2]);

   // then, dot the vorticity with the velocity gradient
   for (i=0;i<nx;i++) {
      for (j=0;j<ny;j++) {
         for (k=0;k<nz;k++) {
            wcgu[0] = wx[i][j][k]*grad[0][0][i][j][k] +
                      wy[i][j][k]*grad[1][0][i][j][k] +
                      wz[i][j][k]*grad[2][0][i][j][k];
            wcgu[1] = wx[i][j][k]*grad[0][1][i][j][k] +
                      wy[i][j][k]*grad[1][1][i][j][k] +
                      wz[i][j][k]*grad[2][1][i][j][k];
            wcgu[2] = wx[i][j][k]*grad[0][2][i][j][k] +
                      wy[i][j][k]*grad[1][2][i][j][k] +
                      wz[i][j][k]*grad[2][2][i][j][k];
            wx[i][j][k] += dt*wcgu[0];
            wy[i][j][k] += dt*wcgu[1];
            wz[i][j][k] += dt*wcgu[2];
         }
      }
      // fprintf(stderr,"%g %g %g\n",wy[i][32][0],wy[i][32][1],wy[i][32][2]);
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
   hxi = 0.5*(float)(nxm1)/SCALE;
   hyi = -0.5*(float)(nxm1);
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
 * create vorticity at solid boundaries - 3D
 *
 * solve 
 */
int create_boundary_vorticity_3d (int nx,int ny,int nz,
      int xbdry,int ybdry,int zbdry,
      float ***u,float ***v,float ***w,
      float ***wx,float ***wy,float ***wz,float dt) {

   int i,j,k;
   int nxm1 = nx-1;
   int nym1 = ny-1;
   int nzm1 = nz-1;
   int nxm2 = nx-2;
   int nym2 = ny-2;
   int nzm2 = nz-2;
   int nxm3 = nx-3;
   int nym3 = ny-3;
   int nzm3 = nz-3;
   float hxi;
   float mult = 1.0*nxm1*dt;	// 1.0 seems to be correct

   // new way, from CottetVM - OOOOooooh, I was using slip walls!
   if (TRUE) {
   hxi = dt*0.5*(float)(nxm1);
   if (xbdry == WALL) {
      for (j=0; j<ny; j++) {
      for (k=0; k<nz; k++) {
         wy[0][j][k] += (w[2][j][k]-4.*w[1][j][k])*hxi;
         wz[0][j][k] += (-1.*v[2][j][k]+4.*v[1][j][k])*hxi;
         wy[nxm1][j][k] += (-1.*w[nxm3][j][k]+4.*w[nxm2][j][k])*hxi;
         wz[nxm1][j][k] += (v[nxm3][j][k]-4.*v[nxm2][j][k])*hxi;
      }
      }
   }
   if (ybdry == WALL) {
      for (i=0; i<nx; i++) {
      for (k=0; k<nz; k++) {
         wx[i][0][k] += (-1.*w[i][2][k]+4.*w[i][1][k])*hxi;
         wz[i][0][k] += (u[i][2][k]-4.*u[i][1][k])*hxi;
         wx[i][nym1][k] += (w[i][nym3][k]-4.*w[i][nym2][k])*hxi;
         wz[i][nym1][k] += (-1.*u[i][nym3][k]+4.*u[i][nym2][k])*hxi;
         // this is the driven cavity!
         // vort[i][nym1][k] += (u[i][nym3][k]-4.*u[i][nym2][k]-3.)*hxi;
      }
      }
   }
   if (zbdry == WALL) {
      for (i=0; i<nx; i++) {
      for (j=0; j<ny; j++) {
         wx[i][j][0] += (v[i][j][2]-4.*v[i][j][1])*hxi;
         wy[i][j][0] += (-1.*u[i][j][2]+4.*u[i][j][1])*hxi;
         wx[i][j][nzm1] += (-1.*v[i][j][nzm3]+4.*v[i][j][nzm2])*hxi;
         wy[i][j][nzm1] += (u[i][j][nzm3]-4.*u[i][j][nzm2])*hxi;
      }
      }
   }
   }

   // my new way (2006) just replace boundary vort with calculated vals
   //   this assumed slip walls, and this routine should be run *before*
   //   the vicmoc advection? (or at least before diffusion)
   if (FALSE) {
   hxi = dt*0.5*(float)(nxm1);
   if (xbdry == WALL) {
     for (j=1; j<nym1; j++) {
     for (k=1; k<nzm1; k++) {
       wx[0][j][k] = (w[0][j+1][k]-w[0][j-1][k])*hxi -
                     (v[0][j][k+1]-v[0][j][k-1])*hxi;
       wx[nxm1][j][k] = (w[nxm1][j+1][k]-w[nxm1][j-1][k])*hxi -
                        (v[nxm1][j][k+1]-v[nxm1][j][k-1])*hxi;
       wy[0][j][k] = (u[0][j][k+1]-u[0][j][k-1])*hxi -
                     (-w[2][j][k]+4.*w[1][j][k]-3.*w[0][j][k])*hxi;
       wy[nxm1][j][k] = (u[nxm1][j][k+1]-u[nxm1][j][k-1])*hxi -
                        (w[nxm1-2][j][k]-4.*w[nxm1-1][j][k]+3.*w[nxm1][j][k])*hxi;
       wz[0][j][k] = (-v[2][j][k]+4.*v[1][j][k]-3.*v[0][j][k])*hxi -
                     (u[0][j+1][k]-u[0][j-1][k])*hxi;
       wz[nxm1][j][k] = (v[nxm1-2][j][k]-4.*v[nxm1-1][j][k]+3.*v[nxm1][j][k])*hxi -
                        (u[nxm1][j+1][k]-u[nxm1][j-1][k])*hxi;
     }
     }
   }
   if (ybdry == WALL) {
     for (i=1; i<nxm1; i++) {
     for (k=1; k<nzm1; k++) {
       wx[i][0][k] = (-w[i][2][k]+4.*w[i][1][k]-3.*w[i][0][k])*hxi -
                     (v[i][0][k+1]-v[i][0][k-1])*hxi;
       wx[i][nym1][k] = (w[i][nym1-2][k]-4.*w[i][nym1-1][k]+3.*w[i][nym1][k])*hxi -
                        (v[i][nym1][k+1]-v[i][nym1][k-1])*hxi;
       wy[i][0][k] = (u[i][0][k+1]-u[i][0][k-1])*hxi -
                     (w[i+1][0][k]-w[i-1][0][k])*hxi;
       wy[i][nym1][k] = (u[i][nym1][k+1]-u[i][nym1][k-1])*hxi -
                        (w[i+1][nym1][k]-w[i-1][nym1][k])*hxi;
       wz[i][0][k] = (v[i+1][0][k]-v[i-1][0][k])*hxi -
                     (-u[i][2][k]+4.*u[i][1][k]-3.*u[i][0][k])*hxi;
       wz[i][nym1][k] = (v[i+1][nym1][k]-v[i-1][nym1][k])*hxi -
                        (u[i][nym1-2][k]-4.*u[i][nym1-1][k]+3.*u[i][nym1][k])*hxi;
     }
     }
   }
   if (zbdry == WALL) {
     for (i=1; i<nxm1; i++) {
     for (j=1; j<nym1; j++) {
       wx[i][j][0] = (w[i][j+1][0]-w[i][j-1][0])*hxi -
                     (-v[i][j][2]+4.*v[i][j][1]-3.*v[i][j][0])*hxi;
       wx[i][j][nzm1] = (w[i][j+1][nzm1]-w[i][j-1][nzm1])*hxi -
                        (v[i][j][nzm1-2]-4.*v[i][j][nzm1-1]+3.*v[i][j][nzm1])*hxi;
       wy[i][j][0] = (-u[i][j][2]+4.*u[i][j][1]-3.*u[i][j][0])*hxi -
                     (w[i+1][j][0]-w[i-1][j][0])*hxi;
       wy[i][j][nzm1] = (u[i][j][nzm1-2]-4.*u[i][j][nzm1-1]+3.*u[i][j][nzm1])*hxi -
                        (w[i+1][j][nzm1]-w[i-1][j][nzm1])*hxi;
       wz[i][j][0] = (v[i+1][j][0]-v[i-1][j][0])*hxi -
                     (u[i][j+1][0]-u[i][j-1][0])*hxi;
       wz[i][j][nzm1] = (v[i+1][j][nzm1]-v[i-1][j][nzm1])*hxi -
                        (u[i][j+1][nzm1]-u[i][j-1][nzm1])*hxi;
     }
     }
   }
   }


   return(0);
}


/* 
 * diffuse scalar field according to laplacian of vorticity, return
 * to a different scalar field
 *
 * where mask is 1.0, do not diffuse, where mask is 0.0, fully diffuse
 *
 * diffus is the scalar diffusivity, like nu for vorticity: 1/Re
 */
float diffuse_scalar_2d (int nx,int ny,int xbdry,int ybdry,float **in,float **out,
      float diffus,int use_MASK,float **mask,float dt) {

   int i,j;
   int nxm1 = nx-1;
   int nym1 = ny-1;
   int istep,numsteps;
   float multx,multy;
   float mult,omm;

   // float mult = 0.001;	// 0.25 is maximum allowed
   // float mult = 0.125;	// 0.25 is maximum allowed
   // float mult = 0.05;	// the new 2nd deriv allows only 0.125 now
   multx = dt*diffus*pow(nxm1/SCALE,2);
   multy = dt*diffus*pow(nxm1,2);

   // if mult>0.125 (stability requirement) then run multiple substeps
   //   each with smaller dt
   if (multx > multy) mult = multx;
   else mult = multy;
   numsteps = 1 + (int)(mult/0.125);
   multx = multx/(float)(numsteps);
   multy = multy/(float)(numsteps);

   for (istep=0; istep<numsteps; istep++) {

   // do the middle part of the field, regardless of the BCs
   //if (use_MASK) {
   if (FALSE) {
     for (i=1;i<nxm1;i++) {
       for (j=1;j<nym1;j++) {
         omm = 1.-mask[i][j];
         out[i][j] = in[i][j] + omm*multx*(in[i+1][j]+in[i-1][j]-2.*in[i][j]) + omm*multy*(in[i][j+1]+in[i][j-1]-2.*in[i][j]);
       }
     }
   } else {
     for (i=1;i<nxm1;i++) {
       for (j=1;j<nym1;j++) {
         out[i][j] = in[i][j] + multx*(in[i+1][j]+in[i-1][j]-2.*in[i][j]) + multy*(in[i][j+1]+in[i][j-1]-2.*in[i][j]);
       }
     }
   }

   // then solve for the wall and periodic bdry
   if (xbdry == WALL || xbdry == OPEN) {
      // for (j=0;j<ny;j++) {
      for (j=1;j<nym1;j++) {
         // out[0][j] = 0.0;
         //out[0][j] = in[0][j] + mult*(in[0][j+1]+in[0][j-1]-in[2][j]+6*in[1][j]-7*in[0][j]);
         out[0][j] = in[0][j] + multx*(-in[2][j]+6.*in[1][j]-5.*in[0][j]) + multy*(in[0][j+1]+in[0][j-1]-2.*in[0][j]);
         // out[nxm1][j] = 0.0;
         //out[nxm1][j] = in[nxm1][j] + mult*(in[nxm1][j+1]+in[nxm1][j-1]-in[nx-3][j]+6*in[nx-2][j]-7*in[nxm1][j]);
         out[nxm1][j] = in[nxm1][j] + multx*(-in[nx-3][j]+6.*in[nx-2][j]-5.*in[nxm1][j]) + multy*(in[nxm1][j+1]+in[nxm1][j-1]-2.*in[nxm1][j]);
      }
   } else if (xbdry == PERIODIC) {
      for (j=1;j<nym1;j++) {
         //out[0][j] = in[0][j] + mult*(in[1][j]+in[nx-2][j]+in[0][j+1]+in[0][j-1]-4.0*in[0][j]);
         out[0][j] = in[0][j] + multx*(in[1][j]+in[nx-2][j]-2.*in[0][j]) + multy*(in[0][j+1]+in[0][j-1]-2.*in[0][j]);
         out[nxm1][j] = out[0][j];
      }
   }
   if (ybdry == WALL || ybdry == OPEN) {
      // for (i=0;i<nx;i++) {
      for (i=1;i<nxm1;i++) {
         // out[i][0] = 0.0;
         //out[i][0] = in[i][0] + mult*(in[i+1][0]+in[i-1][0]-in[i][2]+6*in[i][1]-7*in[i][0]);
         out[i][0] = in[i][0] + multx*(in[i+1][0]+in[i-1][0]-2.*in[i][0]) + multy*(-in[i][2]+6.*in[i][1]-5.*in[i][0]);
         // out[i][nym1] = 0.0;
         //out[i][nym1] = in[i][nym1] + mult*(in[i+1][nym1]+in[i-1][nym1]-in[i][ny-3]+6*in[i][ny-2]-7*in[i][nym1]);
         out[i][nym1] = in[i][nym1] + multx*(in[i+1][nym1]+in[i-1][nym1]-2.*in[i][nym1]) + multy*(-in[i][ny-3]+6.*in[i][ny-2]-5.*in[i][nym1]);
      }
   } else if (ybdry == PERIODIC) {
      for (i=1;i<nxm1;i++) {
         //out[i][0] = in[i][0] + mult*(in[i+1][0]+in[i-1][0]+in[i][1]+in[i][ny-2]-4.0*in[i][0]);
         out[i][0] = in[i][0] + multx*(in[i+1][0]+in[i-1][0]-2.*in[i][0]) + multy*(in[i][1]+in[i][ny-2]-2.*in[i][0]);
         out[i][nym1] = out[i][0];
      }
   }
   if ((xbdry == WALL && ybdry == WALL) || (xbdry == OPEN && ybdry == OPEN)) {
      // corner 2nd derivatives for all 4
      //out[0][0] = in[0][0] + mult*(-in[2][0]-in[0][2]+6*in[1][0]+6*in[0][1]-10*in[0][0]);
      out[0][0] = in[0][0] + multx*(-in[2][0]+6.*in[1][0]-5.*in[0][0]) + multy*(-in[0][2]+6.*in[0][1]-5.*in[0][0]);
      //out[0][nym1] = in[0][nym1] + mult*(-in[2][nym1]-in[0][ny-3]+6*in[1][nym1]+6*in[0][ny-2]-10*in[0][nym1]);
      out[0][nym1] = in[0][nym1] + multx*(-in[2][nym1]+6.*in[1][nym1]-5.*in[0][nym1]) + multy*(-in[0][ny-3]+6.*in[0][ny-2]-5.*in[0][nym1]);
      //out[nxm1][0] = in[nxm1][0] + mult*(-in[nx-3][0]-in[nxm1][2]+6*in[nx-2][0]+6*in[nxm1][1]-10*in[nxm1][0]);
      out[nxm1][0] = in[nxm1][0] + multx*(-in[nx-3][0]+6.*in[nx-2][0]-5.*in[nxm1][0]) + multy*(-in[nxm1][2]+6.*in[nxm1][1]-5.*in[nxm1][0]);
      //out[nxm1][nym1] = in[nxm1][nym1] + mult*(-in[nx-3][nym1]-in[nxm1][ny-3]+6*in[nx-2][nym1]+6*in[nxm1][ny-2]-10*in[nxm1][nym1]);
      out[nxm1][nym1] = in[nxm1][nym1] + multx*(-in[nx-3][nym1]+6.*in[nx-2][nym1]-5.*in[nxm1][nym1]) + multy*(-in[nxm1][ny-3]+6.*in[nxm1][ny-2]-5.*in[nxm1][nym1]);
   } else if (xbdry == WALL && ybdry == PERIODIC) {
      // t-kernel 2nd derivatives for two, then copy
      //out[0][0] = in[0][0] + mult*(in[0][1]+in[0][ny-2]-in[2][0]+6*in[1][0]-7*in[0][0]);
      out[0][0] = in[0][0] + multx*(-in[2][0]+6.*in[1][0]-5.*in[0][0]) + multy*(in[0][1]+in[0][ny-2]-2.*in[0][0]);
      //out[nxm1][0] = in[nxm1][0] + mult*(in[nxm1][1]+in[nxm1][ny-2]-in[nx-3][0]+6*in[nx-2][0]-7*in[nxm1][0]);
      out[nxm1][0] = in[nxm1][0] + multx*(-in[nx-3][0]+6.*in[nx-2][0]-5.*in[nxm1][0]) + multy*(in[nxm1][1]+in[nxm1][ny-2]-2.*in[nxm1][0]);
      out[0][nym1] = out[0][0];
      out[nxm1][nym1] = out[nxm1][0];
   } else if (xbdry == PERIODIC && ybdry == WALL) {
      // t-kernel 2nd derivatives for two, then copy
      //out[0][0] = in[0][0] + mult*(in[1][0]+in[nx-2][0]-in[0][2]+6*in[0][1]-7*in[0][0]);
      out[0][0] = in[0][0] + multx*(in[1][0]+in[nx-2][0]-2.*in[0][0]) + multy*(-in[0][2]+6.*in[0][1]-5.*in[0][0]);
      //out[0][nym1] = in[0][nym1] + mult*(in[1][nym1]+in[nx-2][nym1]-in[0][ny-3]+6*in[0][ny-2]-7*in[0][nym1]);
      out[0][nym1] = in[0][nym1] + multx*(in[1][nym1]+in[nx-2][nym1]-2.*in[0][nym1]) + multy*(-in[0][ny-3]+6.*in[0][ny-2]-5.*in[0][nym1]);
      out[nxm1][0] = out[0][0];
      out[nxm1][nym1] = out[0][nym1];
   } else if (xbdry == PERIODIC && ybdry == PERIODIC) {
      //out[0][0] = in[0][0] + mult*(in[0][1]+in[0][ny-2]+in[1][0]+in[nx-2][0]-4*in[0][0]);
      out[0][0] = in[0][0] + multx*(in[1][0]+in[nx-2][0]-2.*in[0][0]) + multy*(in[0][1]+in[0][ny-2]-2.*in[0][0]);
      out[0][nym1] = out[0][0];
      out[nxm1][0] = out[0][0];
      out[nxm1][nym1] = out[nxm1][0];
   } else {
      fprintf(stderr,"ERROR (diffuse_scalar_2d): should not get here\n");
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

   return(mult*numsteps);
}


/* 
 * diffuse scalar field according to laplacian of vorticity, return
 * to a different scalar field - 3D version
 *
 * diffus is the scalar diffusivity, like nu for vorticity: 1/Re
 */
float diffuse_scalar_3d(int nx,int ny,int nz,int xbdry,int ybdry,int zbdry,
   float ***in,float ***out,float diffus,float dt) {

   int i,j,k;
   int nxm1 = nx-1;
   int nxm2 = nx-2;
   int nym1 = ny-1;
   int nym2 = ny-2;
   int nzm1 = nz-1;
   int nzm2 = nz-2;
   int xs,xf,ys,yf,zs,zf;
   int istep,numsteps;
   float mult;

   // float mult = 0.001;	// 0.25 is maximum allowed
   // float mult = 0.125;	// 0.25 is maximum allowed
   // float mult = 0.05;	// the new 2nd deriv allows only 0.125 now
   mult = dt*diffus*nxm1*nxm1;

   // if mult>0.125 (stability requirement) then run multiple substeps
   //   each with smaller dt
   numsteps = 1 + (int)(mult/0.125);
   mult = mult/(float)(numsteps);

   for (istep=0; istep<numsteps; istep++) {
   // fprintf(stderr,"got here\n"); fflush(stderr);

   // do the middle part of the field, regardless of the BCs
   for (i=1;i<nxm1;i++) {
      for (j=1;j<nym1;j++) {
         for (k=1;k<nzm1;k++) {
            out[i][j][k] = in[i][j][k] + mult*(in[i+1][j][k]+in[i-1][j][k]+in[i][j+1][k]+in[i][j-1][k]+in[i][j][k+1]+in[i][j][k-1]-6.0*in[i][j][k]);
         }
      }
   }
      //printf("  from %g to %g\n",in[8][8][8],out[8][8][8]);
      //printf("    uses %g %g %g %g %g %g\n",in[7][8][8],in[9][8][8],in[8][7][8],in[8][9][8],in[8][8][7],in[8][8][9]);

   // now, do the boundaries, faces first (one bdry variable)
   if (xbdry == WALL) {
      for (j=1;j<nym1;j++) {
      for (k=1;k<nzm1;k++) {
         out[0][j][k] = in[0][j][k] + mult*(
            +in[0][j+1][k]+in[0][j-1][k]
            +in[0][j][k+1]+in[0][j][k-1]
            -in[3][j][k]+4*in[2][j][k]-5*in[1][j][k]
            -2*in[0][j][k]);
            //-in[2][j][k]+6*in[1][j][k]-9*in[0][j][k]);
         out[nxm1][j][k] = in[nxm1][j][k] + mult*(
            +in[nxm1][j+1][k]+in[nxm1][j-1][k]
            +in[nxm1][j][k+1]+in[nxm1][j][k-1]
            -in[nxm1-3][j][k]+4*in[nxm1-2][j][k]-5*in[nxm1-1][j][k]
            -2*in[nxm1][j][k]);
            //-in[nx-3][j][k]+6*in[nx-2][j][k]-9*in[nxm1][j][k]);
      }
      }
   } else if (xbdry == PERIODIC) {
      for (j=1;j<nym1;j++) {
      for (k=1;k<nzm1;k++) {
         out[0][j][k] = in[0][j][k] + mult*(
            +in[0][j+1][k]+in[0][j-1][k]
            +in[0][j][k+1]+in[0][j][k-1]
            +in[1][j][k]+in[nxm1-1][j][k]-6*in[0][j][k]);
            //+in[1][j][k]+in[nx-2][j][k]-6*in[0][j][k]);
         out[nxm1][j][k] = out[0][j][k];
      }
      }
   }
   if (ybdry == WALL) {
      for (i=1;i<nxm1;i++) {
      for (k=1;k<nzm1;k++) {
         out[i][0][k] = in[i][0][k] + mult*(
            +in[i+1][0][k]+in[i-1][0][k]
            +in[i][0][k+1]+in[i][0][k-1]
            -in[i][3][k]+4*in[i][2][k]-5*in[i][1][k]
            -2*in[i][0][k]);
         out[i][nym1][k] = in[i][nym1][k] + mult*(
            +in[i+1][nym1][k]+in[i-1][nym1][k]
            +in[i][nym1][k+1]+in[i][nym1][k-1]
            -in[i][nym1-3][k]+4*in[i][nym1-2][k]-5*in[i][nym1-1][k]
            -2*in[i][nym1][k]);
            //-in[i][ny-3][k]+6*in[i][ny-2][k]-9*in[i][nym1][k]);
      }
      }
      //printf("  from %g to %g\n",in[8][0][8],out[8][0][8]);
      //printf("    uses %g %g %g %g %g %g\n",in[7][0][8],in[9][0][8],in[8][0][8],in[8][1][8],in[8][2][8],in[8][3][8],in[8][0][7],in[8][0][9]);
   } else if (ybdry == PERIODIC) {
      for (i=1;i<nxm1;i++) {
      for (k=1;k<nzm1;k++) {
         out[i][0][k] = in[i][0][k] + mult*(
            +in[i+1][0][k]+in[i-1][0][k]
            +in[i][0][k+1]+in[i][0][k-1]
            +in[i][1][k]+in[i][nym1-1][k]-6*in[i][0][k]);
            //+in[i][1][k]+in[i][ny-2][k]-6*in[i][0][k]);
         out[i][nym1][k] = out[i][0][k];
      }
      }
   }
   if (zbdry == WALL) {
      for (i=1;i<nxm1;i++) {
      for (j=1;j<nym1;j++) {
         out[i][j][0] = in[i][j][0] + mult*(
            +in[i+1][j][0]+in[i-1][j][0]
            +in[i][j+1][0]+in[i][j-1][0]
            -in[i][j][3]+4*in[i][j][2]-5*in[i][j][1]
            -2*in[i][j][0]);
            //-in[i][j][2]+6*in[i][j][1]-9*in[i][j][0]);
         out[i][j][nzm1] = in[i][j][nzm1] + mult*(
            +in[i+1][j][nzm1]+in[i-1][j][nzm1]
            +in[i][j+1][nzm1]+in[i][j-1][nzm1]
            -in[i][j][nzm1-3]+4*in[i][j][nzm1-2]-5*in[i][j][nzm1-1]
            -2*in[i][j][nzm1]);
            //-in[i][j][nz-3]+6*in[i][j][nz-2]-9*in[i][j][nzm1]);
      }
      }
   } else if (zbdry == PERIODIC) {
      for (i=1;i<nxm1;i++) {
      for (j=1;j<nym1;j++) {
         out[i][j][0] = in[i][j][0] + mult*(
            +in[i+1][j][0]+in[i-1][j][0]
            +in[i][j+1][0]+in[i][j-1][0]
            +in[i][j][1]+in[i][j][nzm1-1]-6*in[i][j][0]);
            //+in[i][j][1]+in[i][j][nz-2]-6*in[i][j][0]);
         out[i][j][nzm1] = out[i][j][0];
      }
      }
   }

   // now, the 12 edges! ouch; first the edges in which x varies
   if (ybdry == WALL && zbdry == WALL) {
      for (i=1;i<nxm1;i++) {
         out[i][0][0] = in[i][0][0] + mult*(
            +in[i+1][0][0]+in[i-1][0][0]
            -in[i][2][0]+6*in[i][1][0]
            -in[i][0][2]+6*in[i][0][1]
            -12*in[i][0][0]);
         out[i][0][nzm1] = in[i][0][nzm1] + mult*(
            +in[i+1][0][nzm1]+in[i-1][0][nzm1]
            -in[i][2][nzm1]+6*in[i][1][nzm1]
            -in[i][0][nz-3]+6*in[i][0][nz-2]
            -12*in[i][0][nzm1]);
         out[i][nym1][nzm1] = in[i][nym1][nzm1] + mult*(
            +in[i+1][nym1][nzm1]+in[i-1][nym1][nzm1]
            -in[i][ny-3][nzm1]+6*in[i][ny-2][nzm1]
            -in[i][nym1][nz-3]+6*in[i][nym1][nz-2]
            -12*in[i][nym1][nzm1]);
         out[i][nym1][0] = in[i][nym1][0] + mult*(
            +in[i+1][nym1][0]+in[i-1][nym1][0]
            -in[i][ny-3][0]+6*in[i][ny-2][0]
            -in[i][nym1][2]+6*in[i][nym1][1]
            -12*in[i][nym1][0]);
      }
   } else if (ybdry == PERIODIC && zbdry == WALL) {
      fprintf(stderr,"ERROR (diffuse_scalar_3d): not programmed\n");
      exit(0);
   } else if (ybdry == WALL && zbdry == PERIODIC) {
      fprintf(stderr,"ERROR (diffuse_scalar_3d): not programmed\n");
      exit(0);
   } else if (ybdry == PERIODIC && zbdry == PERIODIC) {
      //fprintf(stderr,"ERROR (diffuse_scalar_3d): not programmed\n");
      //exit(0);
      for (i=1;i<nxm1;i++) {
         out[i][0][0] = in[i][0][0] + mult*(
            +in[i+1][0][0]+in[i-1][0][0]
            +in[i][1][0]+in[i][nym1-1][0]
            +in[i][0][1]+in[i][0][nzm1-1]
            -6*in[i][0][0]);
         out[i][0][nzm1] = out[i][0][0];
         out[i][nym1][nzm1] = out[i][0][0];
         out[i][nym1][0] = out[i][0][0];
      }
   }
   if (xbdry == WALL && zbdry == WALL) {
      for (j=1;j<nym1;j++) {
         out[0][j][0] = in[0][j][0] + mult*(
            +in[0][j+1][0]+in[0][j-1][0]
            -in[2][j][0]+6*in[1][j][0]
            -in[0][j][2]+6*in[0][j][1]
            -12*in[0][j][0]);
         out[0][j][nzm1] = in[0][j][nzm1] + mult*(
            +in[0][j+1][nzm1]+in[0][j-1][nzm1]
            -in[2][j][nzm1]+6*in[1][j][nzm1]
            -in[0][j][nz-3]+6*in[0][j][nz-2]
            -12*in[0][j][nzm1]);
         out[nxm1][j][nzm1] = in[nxm1][j][nzm1] + mult*(
            +in[nxm1][j+1][nzm1]+in[nxm1][j-1][nzm1]
            -in[nx-3][j][nzm1]+6*in[nx-2][j][nzm1]
            -in[nxm1][j][nz-3]+6*in[nxm1][j][nz-2]
            -12*in[nxm1][j][nzm1]);
         out[nxm1][j][0] = in[nxm1][j][0] + mult*(
            +in[nxm1][j+1][0]+in[nxm1][j-1][0]
            -in[nx-3][j][0]+6*in[nx-2][j][0]
            -in[nxm1][j][2]+6*in[nxm1][j][1]
            -12*in[nxm1][j][0]);
      }
   } else if (ybdry == PERIODIC && zbdry == WALL) {
      fprintf(stderr,"ERROR (diffuse_scalar_3d): not programmed\n");
      exit(0);
   } else if (ybdry == WALL && zbdry == PERIODIC) {
      fprintf(stderr,"ERROR (diffuse_scalar_3d): not programmed\n");
      exit(0);
   } else if (ybdry == PERIODIC && zbdry == PERIODIC) {
      //fprintf(stderr,"ERROR (diffuse_scalar_3d): not programmed\n");
      //exit(0);
      for (j=1;j<nym1;j++) {
         out[0][j][0] = in[0][j][0] + mult*(
            +in[0][j+1][0]+in[0][j-1][0]
            +in[1][j][0]+in[nxm1-1][j][0]
            +in[0][j][1]+in[0][j][nzm1-1]
            -6*in[0][j][0]);
         out[0][j][nzm1] = out[0][j][0];
         out[nxm1][j][nzm1] = out[0][j][0];
         out[nxm1][j][0] = out[0][j][0];
      }
   }
   if (xbdry == WALL && ybdry == WALL) {
      for (k=1;k<nzm1;k++) {
         out[0][0][k] = in[0][0][k] + mult*(
            +in[0][0][k+1]+in[0][0][k-1]
            -in[2][0][k]+6*in[1][0][k]
            -in[0][2][k]+6*in[0][1][k]
            -12*in[0][0][k]);
         out[0][nym1][k] = in[0][nym1][k] + mult*(
            +in[0][nym1][k+1]+in[0][nym1][k-1]
            -in[2][nym1][k]+6*in[1][nym1][k]
            -in[0][ny-3][k]+6*in[0][ny-2][k]
            -12*in[0][nym1][k]);
         out[nxm1][nym1][k] = in[nxm1][nym1][k] + mult*(
            +in[nxm1][nym1][k+1]+in[nxm1][nym1][k-1]
            -in[nx-3][nym1][k]+6*in[nx-2][nym1][k]
            -in[nxm1][ny-3][k]+6*in[nxm1][ny-2][k]
            -12*in[nxm1][nym1][k]);
         out[nxm1][0][k] = in[nxm1][0][k] + mult*(
            +in[nxm1][0][k+1]+in[nxm1][0][k-1]
            -in[nx-3][0][k]+6*in[nx-2][0][k]
            -in[nxm1][2][k]+6*in[nxm1][1][k]
            -12*in[nxm1][0][k]);
      }
   } else if (ybdry == PERIODIC && zbdry == WALL) {
      fprintf(stderr,"ERROR (diffuse_scalar_3d): not programmed\n");
      exit(0);
   } else if (ybdry == WALL && zbdry == PERIODIC) {
      fprintf(stderr,"ERROR (diffuse_scalar_3d): not programmed\n");
      exit(0);
   } else if (ybdry == PERIODIC && zbdry == PERIODIC) {
      //fprintf(stderr,"ERROR (diffuse_scalar_3d): not programmed\n");
      //exit(0);
      for (k=1;k<nzm1;k++) {
         out[0][0][k] = in[0][0][k] + mult*(
            +in[0][0][k+1]+in[0][0][k-1]
            +in[1][0][k]+in[nxm1-1][0][k]
            +in[0][1][k]+in[0][nym1-1][k]
            -6*in[0][0][k]);
         out[0][nym1][k] = out[0][0][k];
         out[nxm1][nym1][k] = out[0][0][k];
         out[nxm1][0][k] = out[0][0][k];
      }
   }

   // and, finally, the corners
   if (xbdry == PERIODIC && ybdry == PERIODIC && zbdry == PERIODIC) {
      out[0][0][0] = in[0][0][0] + mult*(
         +in[0][0][1]+in[0][0][nzm1-1]
         +in[1][0][0]+in[nxm1-1][0][0]
         +in[0][1][0]+in[0][nym1-1][0]
         -6*in[0][0][0]);
      out[0][0][nzm1] = out[0][0][0];
      out[0][nym1][0] = out[0][0][0];
      out[0][nym1][nzm1] = out[0][0][0];
      out[nxm1][0][0] = out[0][0][0];
      out[nxm1][0][nzm1] = out[0][0][0];
      out[nxm1][nym1][0] = out[0][0][0];
      out[nxm1][nym1][nzm1] = out[0][0][0];
   }

   // and, for subsequent steps to act on the right variables, swap
   for (i=0;i<nx;i++) {
      for (j=0;j<ny;j++) {
         for (k=0;k<nz;k++) {
            in[i][j][k] = out[i][j][k];
         }
      }
   }

   if (istep > 1) fprintf(stdout,"."); fflush(stdout);
   }	// end loop over numsteps
   if (numsteps > 1) fprintf(stdout,"\n"); fflush(stdout);

   return(mult*numsteps);
}


/*
 * find_vels_2d takes a vorticity field, and solves for the velocity field
 *
 * mask is a flow mask, where 1.0 means solid and 0.0 means open
 *
 * if driven cavity problem, force the lid to vel=1
 */
int find_vels_2d (int step,int isStam,int nx,int ny,int xbdry,int ybdry,
      float *freestream,float **u,float **v,float **vort,int use_MASK,float **mask) {

   int use_multigrid = TRUE;
   int driven_cavity = FALSE;
   int i,j;
   // these are for hwscrt
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
   // these are for mud2sp
   int iparm[17];
   float fparm[6];
   char bndyc[10];
   char cofx[10];
   char cofy[10];
   static float **rhs;
   int mgopt[4];
   static int must_initialize = TRUE;
   // these are for both
   int ierr;
   //static int iworksize = 7000000;	// OK for 768x512
   //static int iworksize = 4000000;	// OK for up to 1024
   //static int iworksize = 16000000;	// OK for 2048
   //static int iworksize = 61655894;	// OK for 4096
   //static int iworksize = 256000000;	// OK for 8192
   //static int iworksize = 384826142;	// OK for 10240
   //static int iworksize = 556000000;	// OK for 12288
   int iworksize = 3.67*nx*ny+100*nx;	// OK for all?
   static float *work;	// increasing this does not fix hwscrt's problem
   float xs = 0.0;
   float ys = 0.0;
   float xf = 1.0;
   //float yf = ((float)(ny-1)/(float)(nx-1));
   //float yf = 1.0;
   float yf = SCALE;
   static float **psi;
   static int allocated_arrays = FALSE;
   static float *left,*right,*top,*bottom;
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
/*
      if (xbdry == WALL) {
         left = (float*)malloc(ny*sizeof(float));
         right = (float*)malloc(ny*sizeof(float));
      }
      if (ybdry == WALL) {
         top = (float*)malloc(nx*sizeof(float));
         bottom = (float*)malloc(nx*sizeof(float));
      }
*/

      // and allocate rhs, if using multigrid solver
      if (use_multigrid)
         rhs = allocate_2d_array_f(nx,ny);

      if (use_MASK) {
         cu = allocate_2d_array_f(nx,ny);
         cv = allocate_2d_array_f(nx,ny);
         cw = allocate_2d_array_f(nx,ny);
         cw_old = allocate_2d_array_f(nx,ny);
      }

      allocated_arrays = TRUE;
   }

   // for all wall boundaries, save the vorticity here, reapply
   //    it later in this routine---we don't want to lose it
/*
   if (xbdry == WALL) {
      for (j=0;j<ny;j++) {
         left[j] = vort[0][j];
         vort[0][j] = 0.0;
         right[j] = vort[nx-1][j];
      }
   }
   if (ybdry == WALL) {
      for (i=0;i<nx;i++) {
         top[i] = vort[i][ny-1];
         bottom[i] = vort[i][0];
      }
   }
*/

   // set maximum number of iterations
   if (use_MASK) {
      maxstep = 100;
   } else {
      maxstep = 1;
   }

   // begin iterative solution
   for (istep=0; istep<maxstep; istep++) {


   // split on solution technique

   // Use mudpack multigrid solver ----------------------------------------
   if (use_multigrid) {

      // set the input parameters
      iparm[0] = 1;	// run a regular solution

      // set the boundary conditions
      if (xbdry == PERIODIC) {
        iparm[1] = 0;
        iparm[2] = 0;
      } else if (xbdry == WALL || xbdry == OPEN) {
        iparm[1] = 1;
        iparm[2] = 1;
        //for (j=0;j<ny;j++) {
        //  psi[0][j] = 0.;
        //  psi[nx-1][j] = 0.;
        //}
      }
      if (ybdry == PERIODIC) {
        iparm[3] = 0;
        iparm[4] = 0;
      } else if (ybdry == WALL || ybdry == OPEN) {
        iparm[3] = 1;
        iparm[4] = 1;
        //for (i=0;i<nx;i++) {
        //  psi[i][0] = 0.;
        //  psi[i][ny-1] = 0.;
        //}
      }

      // find the proper integer dimensions
      // do x first
      //fprintf(stdout," nx is %d\n",nx);
      for (i=2; i<11; i++) {
        for (j=1; j<20; j++) {
          if (i*(pow(2,j))+1 == nx) {
            //fprintf(stdout," nx (%d) is 1 + %d * 2^ %d\n",nx,i,j);
            iparm[5] = i;
            iparm[7] = j+1;
            j = 100;
            i = 100;
          }
        }
      }
      iparm[9] = nx;
      if (i == 11) {
        fprintf(stderr,"Try an nx that is divisible by an integer from");
        fprintf(stderr," 2 to 10\n");
        fprintf(stderr,"Quitting.\n");
        exit(0);
      }
      if (j == 20) {
        fprintf(stderr,"nx = %d ?!? Are you sure?\n",nx);
        fprintf(stderr,"Quitting.\n");
        exit(0);
      }

      // do y next
      for (i=2; i<11; i++) {
        for (j=1; j<20; j++) {
          if (i*(pow(2,j))+1 == ny) {
            //fprintf(stdout," ny (%d) is 1 + %d * 2^ %d\n",ny,i,j);
            iparm[6] = i;
            iparm[8] = j+1;
            j = 100;
            i = 100;
          }
        }
      }
      iparm[10] = ny;
      if (i == 10) {
        fprintf(stderr,"Try an ny that is divisible by an integer from");
        fprintf(stderr," 2 to 10\n");
        fprintf(stderr,"Quitting.\n");
        exit(0);
      }
      if (j == 19) {
        fprintf(stderr,"ny = %d ?!? Are you sure?\n",ny);
        fprintf(stderr,"Quitting.\n");
        exit(0);
      }

      // do some debug printing
      //fprintf(stdout,"iparm(10) is %d\n",iparm[9]);
      //fprintf(stdout,"  ixp %d and iex %d\n",iparm[5],iparm[7]);
      //fprintf(stdout,"iparm(11) is %d\n",iparm[10]);
      //fprintf(stdout,"  jyq %d and jey %d\n",iparm[6],iparm[8]);

      iparm[11] = 0;	// no initial guess for psi is provided
      iparm[12] = 100;	// max # multigrid cycles at finest res
      iparm[13] = 0;	// method of relaxation (0=point, 3=line x and y)
      iparm[14] = iworksize;	// size of workspace

      fparm[0] = xs;	// x start
      fparm[1] = xf;	// x end
      fparm[2] = ys;	// y start
      fparm[3] = yf;	// y end
      fparm[4] = 1.e-4;	// tolerance criterion

      mgopt[0] = 0;	// use default multigrid options

      // the boundary condition subroutine names
      // Apparently, the F77 solver calls these as "external"
      // I have no clue how to actually get them to work
      // For now, I have just replaced the calls in mud2sp_full.f with val=1.
      (void) strcpy(cofx,"cofx");
      (void) strcpy(cofy,"cofy");
      // bndyc only used for mixed boundary conditions - test Dirichlet conditions first
      (void) strcpy(bndyc,"bndyc");

      // only do this once
      if (must_initialize) {

         // this flags the routine to discretize the problem and check input
         iparm[0] = 0;
         fprintf(stdout,"Running initialization mud2sp\n");
         mud2sp_(iparm, fparm, work, cofx, cofy,
                 bndyc, rhs[0], psi[0], mgopt, &ierr);

         // catch an error in the setup or input
         if (ierr != 0) {
            fprintf(stderr,"ERROR (mud2sp_): ierr = %d\n",ierr);
            if (ierr == 9) fprintf(stderr,"  iparm[14] too small, > %d\n",iparm[15]);
            fprintf(stderr,"Quitting.\n");
            exit(0);
         }

         // read workspace requirements!
         //iworksize = iparm[15]+100;
         //fprintf(stderr,"work array needs %d\n",iworksize);
         fprintf(stderr,"work array needs %d\n",iparm[15]);
         fprintf(stderr,"  we currently allocate %d\n",iworksize);
         //free_1d_array_f(work);
         //work = allocate_1d_array_f(iworksize);

         // don't do this again
         must_initialize = FALSE;
      }

      // now call the subroutine for real

      // set the rhs values (negative vorticity)
      for (i=0;i<nx;i++) {
         for (j=0;j<ny;j++) {
            rhs[i][j] = -1.0*vort[i][j];
            //psi[i][j] = 0.;
            //if (abs(rhs[i][j]) > 20.) fprintf(stderr,"rhs %d %d too large %g\n",i,j,rhs[i][j]);
         }
      }

      // initialize psi
      for (i=0;i<nx;i++) {
         for (j=0;j<ny;j++) {
            psi[i][j] = 0.;
         }
      }

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

      iparm[0] = 1;	// run a regular solution

      //fprintf(stdout,"Running mud2sp solution\n");
      fprintf(stderr,"    running mud2sp\n"); fflush(stderr);
      mud2sp_(iparm, fparm, work, cofx, cofy,
              bndyc, rhs[0], psi[0], mgopt, &ierr);

      //mud2sp_(int *iparm, float *fparm, float *work, char *cfx, char* cfy,
      //        char *bndyc, float *rhs, float *phi, int *mgopt, int ierror);

      //fprintf(stdout,"mud2sp solved in %d cycles\n",iparm[16]);

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
   //sprintf(outfileroot,"psi_%04d",step);
   //write_output(outfileroot,nx,ny,psi,-0.1,0.2,1);

   // now, find the derivatives of streamfunction to make the velocities!
   find_gradient_of_scalar_2nd_2d (nx,ny,xbdry,ybdry,psi,v,-1.0,u,1.0);

   // if there's a mask, then find the counter-vorticity
   if (use_MASK) {

      // first, find counter-velocity
      for (i=0;i<nx;i++) {
         for (j=0;j<ny;j++) {
            cu[i][j] = -1.0*u[i][j]*mask[i][j];
            cv[i][j] = -1.0*v[i][j]*mask[i][j];
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
            vort[i][j] = (1.-mask[i][j]) * (vort[i][j] + cw[i][j]);
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
         if (maxcorr/maxvort < 1.e-3) {
            istep += maxstep;
         }
         if (isnan(maxcorr/maxvort)) {
           fprintf(stderr,"Quitting.\n");
           exit(1);
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
 * find_vels_3d takes a vorticity field, and solves for the velocity field - in 3D
 */
int find_vels_3d (int step,int nx,int ny,int nz,int xbdry,int ybdry,int zbdry,
   float ***u,float ***v,float ***w,float ***wx,float ***wy,float ***wz) {

   int use_multigrid = FALSE;
   int slip_wall = FALSE;
   int driven_cavity = FALSE;
   int i,j,k;
   int nxm1 = nx-1;
   int nym1 = ny-1;
   int nzm1 = nz-1;
   int nxp = nx;
   int nyp = ny;
   int nzp = nz;
   int bcux,bcvx,bcwx,bcuy,bcvy,bcwy,bcuz,bcvz,bcwz;
   float elmbda = 0.0;
   // these are for mud2sp
   int iparm[22];
   float fparm[8];
   char bndyc[10];
   char cofx[10];
   char cofy[10];
   char cofz[10];
   static float ***rhs;
   int mgopt[4];
   // these are for both
   int ierr;
   static int iworksize = 4000000;
   static float *work;	// increasing this does not fix hwscrt's problem
   static int allocated_arrays = FALSE;
   static int must_initialize = TRUE;
   float pert;
   float xs = 0.0;
   float ys = 0.0;
   float zs = 0.0;
   float xf = 1.0;
   float yf = (float)((float)(ny-1)/(float)(nx-1));
   float zf = (float)((float)(nz-1)/(float)(nx-1));
   static int malloced_yet = FALSE;
   static float **bduxs,**bduxf,**bduys,**bduyf,**bduzs,**bduzf;
   static float **bdvxs,**bdvxf,**bdvys,**bdvyf,**bdvzs,**bdvzf;
   static float **bdwxs,**bdwxf,**bdwys,**bdwyf,**bdwzs,**bdwzf;
   static float **xswx,**xswy,**xswz,**xfwx,**xfwy,**xfwz;
   static float **yswx,**yswy,**yswz,**yfwx,**yfwy,**yfwz;
   static float **zswx,**zswy,**zswz,**zfwx,**zfwy,**zfwz;
   static float ***grad[3][3];
   //static float ***div;
   float div,maxdiv;
   static int set_grad = TRUE;		// FALSE if you need it
   char outfileroot[70];

   // fprintf(stderr,"in find_vels_3d\n"); fflush(stderr);

   // allocate memory for the streamfunction, and the boundary
   //    vorticity arrays
   if (!allocated_arrays) {

      // always allocate work array
      work = allocate_1d_array_f(iworksize);

      // and allocate rhs, if using multigrid solver
      if (use_multigrid)
         rhs = allocate_3d_array_f(nx,ny,nz);

      allocated_arrays = TRUE;
   }

   if (!use_multigrid) {

   // allocate memory for the boundary condition arrays
   if (!malloced_yet) {
      // if (xbdry == WALL) {
         bduxs = allocate_2d_array_f(ny,nz);
         bduxf = allocate_2d_array_f(ny,nz);
      // }
      // if (ybdry == WALL) {
         bduys = allocate_2d_array_f(nx,nz);
         bduyf = allocate_2d_array_f(nx,nz);
      // }
      // if (zbdry == WALL) {
         bduzs = allocate_2d_array_f(nx,ny);
         bduzf = allocate_2d_array_f(nx,ny);
      // }
      bdvxs = allocate_2d_array_f(ny,nz);
      bdvxf = allocate_2d_array_f(ny,nz);
      bdvys = allocate_2d_array_f(nx,nz);
      bdvyf = allocate_2d_array_f(nx,nz);
      bdvzs = allocate_2d_array_f(nx,ny);
      bdvzf = allocate_2d_array_f(nx,ny);
      bdwxs = allocate_2d_array_f(ny,nz);
      bdwxf = allocate_2d_array_f(ny,nz);
      bdwys = allocate_2d_array_f(nx,nz);
      bdwyf = allocate_2d_array_f(nx,nz);
      bdwzs = allocate_2d_array_f(nx,ny);
      bdwzf = allocate_2d_array_f(nx,ny);
      malloced_yet = TRUE;

      xswx = allocate_2d_array_f(ny,nz);
      xswy = allocate_2d_array_f(ny,nz);
      xswz = allocate_2d_array_f(ny,nz);
      xfwx = allocate_2d_array_f(ny,nz);
      xfwy = allocate_2d_array_f(ny,nz);
      xfwz = allocate_2d_array_f(ny,nz);

      yswx = allocate_2d_array_f(nx,nz);
      yswy = allocate_2d_array_f(nx,nz);
      yswz = allocate_2d_array_f(nx,nz);
      yfwx = allocate_2d_array_f(nx,nz);
      yfwy = allocate_2d_array_f(nx,nz);
      yfwz = allocate_2d_array_f(nx,nz);

      zswx = allocate_2d_array_f(nx,ny);
      zswy = allocate_2d_array_f(nx,ny);
      zswz = allocate_2d_array_f(nx,ny);
      zfwx = allocate_2d_array_f(nx,ny);
      zfwy = allocate_2d_array_f(nx,ny);
      zfwz = allocate_2d_array_f(nx,ny);
   }

   // allocate memory for the gradient
   if (!set_grad) {
      grad[0][0] = allocate_3d_array_f(nx,ny,nz);
      grad[0][1] = allocate_3d_array_f(nx,ny,nz);
      grad[0][2] = allocate_3d_array_f(nx,ny,nz);
      grad[1][0] = allocate_3d_array_f(nx,ny,nz);
      grad[1][1] = allocate_3d_array_f(nx,ny,nz);
      grad[1][2] = allocate_3d_array_f(nx,ny,nz);
      grad[2][0] = allocate_3d_array_f(nx,ny,nz);
      grad[2][1] = allocate_3d_array_f(nx,ny,nz);
      grad[2][2] = allocate_3d_array_f(nx,ny,nz);
      //div = allocate_3d_array_f(nx,ny,nz);
      set_grad = TRUE;
   }

   // if a boundary is periodic, make sure that both ends have the same value

   // try this: save the vorticity on the z-min boundary and
   //  zero it before taking the curl
   // THIS WORKS
   // whether slip or no-slip walls, do this:

   if (FALSE) {
      if (xbdry == WALL) {
         for (j=0;j<ny;j++) {
            for (k=0;k<nz;k++) {
               xswx[j][k] = wx[0][j][k];
               xswy[j][k] = wy[0][j][k];
               xswz[j][k] = wz[0][j][k];
               wx[0][j][k] = 0.;
               wy[0][j][k] = 0.;
               wz[0][j][k] = 0.;
               xfwx[j][k] = wx[nxm1][j][k];
               xfwy[j][k] = wy[nxm1][j][k];
               xfwz[j][k] = wz[nxm1][j][k];
               wx[nxm1][j][k] = 0.;
               wy[nxm1][j][k] = 0.;
               wz[nxm1][j][k] = 0.;
            }
         }
      }
      if (ybdry == WALL) {
         for (i=0;i<nx;i++) {
            for (k=0;k<nz;k++) {
               yswx[i][k] = wx[i][0][k];
               yswy[i][k] = wy[i][0][k];
               yswz[i][k] = wz[i][0][k];
               wx[i][0][k] = 0.;
               wy[i][0][k] = 0.;
               wz[i][0][k] = 0.;
               yfwx[i][k] = wx[i][nym1][k];
               yfwy[i][k] = wy[i][nym1][k];
               yfwz[i][k] = wz[i][nym1][k];
               wx[i][nym1][k] = 0.;
               wy[i][nym1][k] = 0.;
               wz[i][nym1][k] = 0.;
            }
         }
      }
      if (zbdry == WALL) {
         for (i=0;i<nx;i++) {
            for (j=0;j<ny;j++) {
               zswx[i][j] = wx[i][j][0];
               zswy[i][j] = wy[i][j][0];
               zswz[i][j] = wz[i][j][0];
               wx[i][j][0] = 0.;
               wy[i][j][0] = 0.;
               wz[i][j][0] = 0.;
               zfwx[i][j] = wx[i][j][nzm1];
               zfwy[i][j] = wy[i][j][nzm1];
               zfwz[i][j] = wz[i][j][nzm1];
               wx[i][j][nzm1] = 0.;
               wy[i][j][nzm1] = 0.;
               wz[i][j][nzm1] = 0.;
            }
         }
      }
   }

   } // end if not multigrid

   // initialize the right-hand-side as the negative curl of the vorticity
   // one way (1st order on boundaries)
   compute_curl_3d(nx,ny,nz,xbdry,ybdry,zbdry,wx,wy,wz,u,v,w);

   // other way (2nd order globally)
   // find_gradient_of_scalar_3d(nx,ny,nz,xbdry,ybdry,zbdry,wx,grad[0]);
   // find_gradient_of_scalar_3d(nx,ny,nz,xbdry,ybdry,zbdry,wy,grad[1]);
   // find_gradient_of_scalar_3d(nx,ny,nz,xbdry,ybdry,zbdry,wz,grad[2]);

   // out of curiosity, what is the vorticity divergence?
   //div = allocate_3d_array_f(nx,ny,nz);
   //maxdiv=0.;
   //for (i=1; i<nxm1; i++)
   //for (j=1; j<nym1; j++)
   //for (k=1; k<nzm1; k++) {
   //  div = 0.5*nxm1*(wx[i+1][j][k]-wx[i-1][j][k]+wy[i][j+1][k]-wy[i][j-1][k]+wz[i][j][k+1]-wz[i][j][k-1]);
   //  if (div>maxdiv) {
   //    maxdiv=div;
       //printf("    %g at %d %d %d\n",div,i,j,k);
       //printf("          %g %g %g %g %g %g\n",wx[i+1][j][k],wx[i-1][j][k],wy[i][j+1][k],wy[i][j-1][k],wz[i][j][k+1],wz[i][j][k-1]);
     //}
     //div[i][j][k] = grad[0][0][i][j][k]+grad[1][1][i][j][k]+grad[2][2][i][j][k];
   //}
   //printf("  max div %g\n",maxdiv);
   // sprintf(outfileroot,"div_vort_%04d",step);
   // write_output_3d(outfileroot,nx,ny,nz,div,-0.01,0.02,2,FALSE);


   // and make sure that its negative
   for (i=0; i<nx; i++)
   for (j=0; j<ny; j++)
   for (k=0; k<nz; k++) {
      u[i][j][k] = -u[i][j][k];
      v[i][j][k] = -v[i][j][k];
      w[i][j][k] = -w[i][j][k];
   }

   // split on solution technique

   // Use mudpack multigrid solver ----------------------------------------
   if (use_multigrid) {

      // set the input parameters
      iparm[0] = 1;	// run a regular solution

      // set the boundary conditions
      if (xbdry == PERIODIC) {
        iparm[1] = 0;
        iparm[2] = 0;
      } else {
        iparm[1] = 1;
        iparm[2] = 1;
        for (j=0;j<ny;j++) {
        for (k=0;k<nz;k++) {
          //psi[0][j] = 0.;
          //psi[nx-1][j] = 0.;
        }
        }
      }
      if (ybdry == PERIODIC) {
        iparm[3] = 0;
        iparm[4] = 0;
      } else {
        iparm[3] = 1;
        iparm[4] = 1;
        for (i=0;i<nx;i++) {
          //psi[i][0] = 0.;
          //psi[i][ny-1] = 0.;
        }
      }
      if (zbdry == PERIODIC) {
        iparm[5] = 0;
        iparm[6] = 0;
      } else {
        iparm[5] = 1;
        iparm[6] = 1;
        for (i=0;i<nx;i++) {
          //psi[i][0] = 0.;
          //psi[i][ny-1] = 0.;
        }
      }

      // find the proper integer dimensions
      // do x first
      for (i=2; i<11; i++) {
        for (j=1; j<20; j++) {
          if (i*(pow(2,j-1))+1 == nx) {
            //fprintf(stdout," nx is %d * %d^2 + 1\n",i,j);
            iparm[7] = i;
            iparm[10] = j;
            j = 100;
            i = 100;
          }
        }
      }
      iparm[13] = nx;
      if (i == 11) {
        fprintf(stderr,"Try an nx that is divisible by an integer from");
        fprintf(stderr," 2 to 10\n");
        fprintf(stderr,"Quitting.\n");
        exit(0);
      }
      if (j == 20) {
        fprintf(stderr,"nx = %d ?!? Are you sure?\n",nx);
        fprintf(stderr,"Quitting.\n");
        exit(0);
      }

      // do y next
      for (i=2; i<11; i++) {
        for (j=1; j<20; j++) {
          if (i*(pow(2,j-1))+1 == ny) {
            //fprintf(stdout," ny is %d * %d^2 + 1\n",i,j);
            iparm[8] = i;
            iparm[11] = j;
            j = 100;
            i = 100;
          }
        }
      }
      iparm[14] = ny;
      if (i == 10) {
        fprintf(stderr,"Try an ny that is divisible by an integer from");
        fprintf(stderr," 2 to 10\n");
        fprintf(stderr,"Quitting.\n");
        exit(0);
      }
      if (j == 19) {
        fprintf(stderr,"ny = %d ?!? Are you sure?\n",ny);
        fprintf(stderr,"Quitting.\n");
        exit(0);
      }

      // do z last
      for (i=2; i<11; i++) {
        for (j=1; j<20; j++) {
          if (i*(pow(2,j-1))+1 == ny) {
            //fprintf(stdout," ny is %d * %d^2 + 1\n",i,j);
            iparm[9] = i;
            iparm[12] = j;
            j = 100;
            i = 100;
          }
        }
      }
      iparm[15] = ny;

      // do some debug printing
      //fprintf(stdout,"iparm(10) is %d\n",iparm[9]);
      //fprintf(stdout,"  ixp %d and iex %d\n",iparm[5],iparm[7]);
      //fprintf(stdout,"iparm(11) is %d\n",iparm[10]);
      //fprintf(stdout,"  jyq %d and jey %d\n",iparm[6],iparm[8]);

      iparm[16] = 0;	// no initial guess for psi is provided
      iparm[17] = 100;	// max # multigrid cycles at finest res
      iparm[18] = 0;	// method of relaxation (0=point, 3=line x and y)
      iparm[19] = iworksize;	// size of workspace

      fparm[0] = 0.0;	// x start
      fparm[1] = 1.0;	// x end
      fparm[2] = 0.0;	// y start
      fparm[3] = 1.0;	// y end
      fparm[4] = 0.0;	// z start
      fparm[5] = 1.0;	// z end
      fparm[6] = 1.e-4;	// tolerance criterion

      mgopt[0] = 0;	// use default multigrid options

      // the boundary condition subroutine names
      // Apparently, the F77 solver calls these as "external"
      // I have no clue how to actually get them to work
      // For now, I have just replaced the calls in mud2sp_full.f with val=1.
      (void) strcpy(cofx,"dummy");
      (void) strcpy(cofy,"dummy");
      (void) strcpy(cofz,"dummy");
      // bndyc only used for mixed boundary conditions - test Dirichlet conditions first
      (void) strcpy(bndyc,"dummy");

      // only do this once
      if (must_initialize) {

         // this flags the routine to discretize the problem and check input
         iparm[0] = 0;
         fprintf(stdout,"Running initialization mud3sp\n");
#ifdef IS_3D
         mud3sp_(iparm, fparm, work, cofx, cofy, cofz,
                 bndyc, rhs[0][0], u[0][0], mgopt, &ierr);
#endif

         // catch an error in the setup or input
         if (ierr != 0) {
            fprintf(stderr,"ERROR (mud3sp_): ierr = %d\n",ierr);
            if (ierr == 9) fprintf(stderr,"  iparm[19] too small, > %d\n",iparm[20]);
            fprintf(stderr,"Quitting.\n");
            exit(0);
         }

         // read workspace requirements!
         //iworksize = iparm[15]+100;
         //fprintf(stderr,"work array needs %d\n",iworksize);
         //free_1d_array_f(work);
         //work = allocate_1d_array_f(iworksize);


         // don't do this again
         must_initialize = FALSE;
      }

      // now call the subroutine for real

//    I'M HERE  ----------------------------------------------------------------------------**

      // set the rhs values (negative vorticity)
      for (i=0;i<nx;i++) {
        for (j=0;j<ny;j++) {
          for (j=0;j<ny;j++) {
            //rhs[i][j][k] = -1.0*vort[i][j];
            //psi[i][j] = 0.;
            //if (abs(rhs[i][j]) > 20.) fprintf(stderr,"rhs %d %d too large %g\n",i,j,rhs[i][j]);
          }
        }
      }

      // initialize psi
      for (i=0;i<nx;i++) {
        for (j=0;j<ny;j++) {
          //psi[i][j] = 0.;
        }
      }

      iparm[0] = 1;	// run a regular solution

      //fprintf(stdout,"Running mud3sp solution\n");
#ifdef IS_3D
      mud3sp_(iparm, fparm, work, cofx, cofy, cofz,
              bndyc, rhs[0][0], u[0][0], mgopt, &ierr);
#endif

      //mud2sp_(int *iparm, float *fparm, float *work, char *cfx, char* cfy,
      //        char *bndyc, float *rhs, float *phi, int *mgopt, int ierror);

      //fprintf(stdout,"mud2sp solved in %d cycles\n",iparm[16]);

      // debug print the resulting streamfunctions
/*
      for (i=0; i<8; i++) {
         for (j=0; j<8; j++) {
            fprintf(stdout," %g",psi[i][j]);
         }
         fprintf(stdout,"\n");
      }
*/

      // catch a runtime error
      if (ierr != 0) {
         fprintf(stderr,"ERROR (mud3sp_): ierr = %d\n",ierr);
         fprintf(stderr,"Quitting.\n");
         exit(0);
      }

      //exit(0);

   // Use fishpak FFT solver ----------------------------------------------
   } else {

   // for all wall boundaries, set zero thru-flow and derivative
   if (slip_wall) {
      // allow slip at walls (LES-like)
      if (xbdry == WALL) {
         for (j=0;j<ny;j++) {
            for (k=0;k<nz;k++) {
               u[0][j][k] = 0.0;
               u[nxm1][j][k] = 0.0;
               bdvxs[j][k] = 0.0;
               bdvxf[j][k] = 0.0;
               bdwxs[j][k] = 0.0;
               bdwxf[j][k] = 0.0;
            }
         }
      }
      if (ybdry == WALL) {
         for (i=0;i<nx;i++) {
            for (k=0;k<nz;k++) {
               v[i][0][k] = 0.0;
               v[i][nym1][k] = 0.0;
               bduys[i][k] = 0.0;
               bduyf[i][k] = 0.0;
               bdwys[i][k] = 0.0;
               bdwyf[i][k] = 0.0;
            }
         }
      }
      if (zbdry == WALL) {
         for (i=0;i<nx;i++) {
            for (j=0;j<ny;j++) {
               w[i][j][0] = 0.0;
               w[i][j][nzm1] = 0.0;
               bduzs[i][j] = 0.0;
               bduzf[i][j] = 0.0;
               bdvzs[i][j] = 0.0;
               bdvzf[i][j] = 0.0;
            }
         }
      }
   } else {
      // no-slip wall (viscous simulation)
      if (xbdry == WALL) {
         for (j=0;j<ny;j++) {
            for (k=0;k<nz;k++) {
               u[0][j][k] = 0.0;
               u[nxm1][j][k] = 0.0;
               // and tangential components should be zero, too
               v[0][j][k] = 0.0;
               v[nxm1][j][k] = 0.0;
               w[0][j][k] = 0.0;
               w[nxm1][j][k] = 0.0;
            }
         }
      }
      if (ybdry == WALL) {
         for (i=0;i<nx;i++) {
            for (k=0;k<nz;k++) {
               v[i][0][k] = 0.0;
               v[i][nym1][k] = 0.0;
               // and tangential components should be zero, too
               u[i][0][k] = 0.0;
               u[i][nym1][k] = 0.0;
               w[i][0][k] = 0.0;
               w[i][nym1][k] = 0.0;
            }
         }
      }
      if (zbdry == WALL) {
         for (i=0;i<nx;i++) {
            for (j=0;j<ny;j++) {
               w[i][j][0] = 0.0;
               w[i][j][nzm1] = 0.0;
               // and tangential components should be zero, too
               u[i][j][0] = 0.0;
               u[i][j][nzm1] = 0.0;
               v[i][j][0] = 0.0;
               v[i][j][nzm1] = 0.0;
            }
         }
      }
   }

   // force the driven cavity using boundary conditions!
   if (driven_cavity) {
      for (i=0;i<nx;i++)
      for (j=0;j<ny;j++) {
         // dummy
      }
   }

   // write out RHS
   // sprintf(outfileroot,"dcvx_%04d",step);
   // write_output_3d(outfileroot,nx,ny,nz,u,-500,1000,2,FALSE);
   // sprintf(outfileroot,"dcvy_%04d",step);
   // write_output_3d(outfileroot,nx,ny,nz,v,-500,1000,2,FALSE);
   // sprintf(outfileroot,"dcvz_%04d",step);
   // write_output_3d(outfileroot,nx,ny,nz,w,-500,1000,2,FALSE);

   // set bc[uvw][xyz] flag boundary conditions
   // for '1' (value specified at boundary), the value is in u/v/w,
   //   and bd[xyz][sf] is not even used (no-slip walls)
   // for '3' (derivative specified at bdry), u/v/w is the rhs, and
   //   bd[xyz][sf] is the derivative (slip wall case)
   if (xbdry == PERIODIC) {
      bcux = 0;
      bcvx = 0;
      bcwx = 0;
   } else if (xbdry == WALL) {
      /* both sides are walls */
      if (slip_wall) {
         bcux = 1;
         bcvx = 3;
         bcwx = 3;
      } else {
         bcux = 1;
         bcvx = 1;
         bcwx = 1;
      }
   } else {
      fprintf(stdout,"ERROR (find_vels_3d): BCs for x-dir invalid\n");
      fprintf(stdout,"exiting.\n");
      exit(0);
   }
   if (ybdry == PERIODIC) {
      bcuy = 0;
      bcvy = 0;
      bcwy = 0;
   } else if (ybdry == WALL) {
      /* both sides are walls */
      if (slip_wall) {
         bcuy = 3;
         bcvy = 1;
         bcwy = 3;
      } else {
         bcuy = 1;
         bcvy = 1;
         bcwy = 1;
      }
   } else {
      fprintf(stdout,"ERROR (find_vels_3d): BCs for y-dir invalid\n");
      fprintf(stdout,"exiting.\n");
      exit(0);
   }
   if (zbdry == PERIODIC) {
      bcuz = 0;
      bcvz = 0;
      bcwz = 0;
   } else if (zbdry == WALL) {
      /* both sides are walls */
      if (slip_wall) {
         bcuz = 3;
         bcvz = 3;
         bcwz = 1;
      } else {
         bcuz = 1;
         bcvz = 1;
         bcwz = 1;
      }
   } else {
      fprintf(stdout,"ERROR (find_vels_3d): BCs for z-dir invalid\n");
      fprintf(stdout,"exiting.\n");
      exit(0);
   }

   // make the call to the external solver, note that x and y are backwards...

   // put in cvx as the 3D array, it returns velocity in its place
   // hw3crt_(&f->z.s[2],&f->z.f[2],&f->z.n[2],&f->bcuz,f->bdxs[0],f->bdxf[0],
   //         &f->z.s[1],&f->z.f[1],&f->z.n[1],&f->bcuy,f->bdys[0],f->bdyf[0],
   //         &f->z.s[0],&f->z.f[0],&f->z.n[0],&f->bcux,f->bdzs[0],f->bdzf[0],
   //         &elmbda,&f->z.np1[2],&f->z.np1[1],f->fi[at]->u[0][0],&pert,&ierr,work);

// must comment these out to get 2D code working
#ifdef IS_3D
   hw3crt_(&zs,&zf,&nzm1,&bcuz,bduxs[0],bduxf[0],
           &ys,&yf,&nym1,&bcuy,bduys[0],bduyf[0],
           &xs,&xf,&nxm1,&bcux,bduzs[0],bduzf[0],
           &elmbda,&nzp,&nyp,u[0][0],&pert,&ierr,work);
   hw3crt_(&zs,&zf,&nzm1,&bcvz,bdvxs[0],bdvxf[0],
           &ys,&yf,&nym1,&bcvy,bdvys[0],bdvyf[0],
           &xs,&xf,&nxm1,&bcvx,bdvzs[0],bdvzf[0],
           &elmbda,&nzp,&nyp,v[0][0],&pert,&ierr,work);
   hw3crt_(&zs,&zf,&nzm1,&bcwz,bdwxs[0],bdwxf[0],
           &ys,&yf,&nym1,&bcwy,bdwys[0],bdwyf[0],
           &xs,&xf,&nxm1,&bcwx,bdwzs[0],bdwzf[0],
           &elmbda,&nzp,&nyp,w[0][0],&pert,&ierr,work);
#endif

   // fprintf(stdout,"  hwscrt used %d locations in work array, err=%d\n",(int)(work[0]),ierr);

   // find vel gradients
   // find_gradient_of_scalar_3d(nx,ny,nz,xbdry,ybdry,zbdry,u,grad[0]);
   // find_gradient_of_scalar_3d(nx,ny,nz,xbdry,ybdry,zbdry,v,grad[1]);
   // find_gradient_of_scalar_3d(nx,ny,nz,xbdry,ybdry,zbdry,w,grad[2]);

   // out of curiosity, what is the vorticity divergence?
   // for (i=0; i<nx; i++)
   // for (j=0; j<ny; j++)
   // for (k=0; k<nz; k++) {
      // div[i][j][k] = grad[0][0][i][j][k]+grad[1][1][i][j][k]+grad[2][2][i][j][k];
   // }
   // sprintf(outfileroot,"div_vel_%04d",step);
   // write_output_3d(outfileroot,nx,ny,nz,div,-1.,2.,2,FALSE);

   // and what is wy?
   // for (i=0; i<nx; i++)
   // for (j=0; j<ny; j++)
   // for (k=0; k<nz; k++) {
      // div[i][j][k] = grad[0][2][i][j][k]-grad[2][0][i][j][k];
   // }
   // sprintf(outfileroot,"wy_after_%04d",step);
   // write_output_3d(outfileroot,nx,ny,nz,div,-10.,20.,2,FALSE);

   // put the vorticity that we removed back
   if (FALSE) {
      if (xbdry == WALL) {
         for (j=0;j<ny;j++) {
            for (k=0;k<nz;k++) {
               wx[0][j][k] = xswx[j][k];
               wy[0][j][k] = xswy[j][k];
               wz[0][j][k] = xswz[j][k];
               wx[nxm1][j][k] = xfwx[j][k];
               wy[nxm1][j][k] = xfwy[j][k];
               wz[nxm1][j][k] = xfwz[j][k];
            }
         }
      }
      if (ybdry == WALL) {
         for (i=0;i<nx;i++) {
            for (k=0;k<nz;k++) {
               wx[i][0][k] = yswx[i][k];
               wy[i][0][k] = yswy[i][k];
               wz[i][0][k] = yswz[i][k];
               wx[i][nym1][k] = yfwx[i][k];
               wy[i][nym1][k] = yfwy[i][k];
               wz[i][nym1][k] = yfwz[i][k];
            }
         }
      }
      if (zbdry == WALL) {
         for (i=0;i<nx;i++) {
            for (j=0;j<ny;j++) {
               wx[i][j][0] = zswx[i][j];
               wy[i][j][0] = zswy[i][j];
               wz[i][j][0] = zswz[i][j];
               wx[i][j][nzm1] = zfwx[i][j];
               wy[i][j][nzm1] = zfwy[i][j];
               wz[i][j][nzm1] = zfwz[i][j];
            }
         }
      }
   }

   } // end if not multigrid

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
   if (nx>ny) vel = allocate_1d_array_f(nx);
   else vel = allocate_1d_array_f(ny);

   // initialize the static variables
   find_biot_savart (0.,0.,nx,ny,-1.0,-1.0,vort,thisvel);

   // first, march across bottom row (j=0, y=0)
   #pragma omp parallel for private(i,thisvel)
   for (i=0;i<nx;i++) {
     find_biot_savart (i*dx,0.,nx,ny,dx,dy,vort,thisvel);
     vel[i] = freestream[1] + thisvel[1];
     //vel[i] = freestream[1] + find_biot_savart_v(i*dx,0.,nx,ny,dx,dy,vort);
     //fprintf(stderr,"fbs %g, fbsv %g, n %d\n",thisvel[1],vel[i]-freestream[1],j);
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
     find_biot_savart (0.,j*dy,nx,ny,dx,dy,vort,thisvel);
     vel[j] = freestream[0] + thisvel[0];
     //vel[j] = freestream[0] + find_biot_savart_u(0.,j*dy,nx,ny,dx,dy,vort);
   }
   for (j=1;j<ny;j++) {
     psi[0][j] = psi[0][j-1] + 0.5*(vel[j]+vel[j-1])*dy;
     //fprintf(stderr,"y=%g  u=%g  psi=%g\n",j*dy,vel[j],psi[0][j]);
   }

   // then march up right side
   #pragma omp parallel for private(j,thisvel)
   for (j=0;j<ny;j++) {
     find_biot_savart (1.,j*dy,nx,ny,dx,dy,vort,thisvel);
     vel[j] = freestream[0] + thisvel[0];
     //vel[j] = freestream[0] + find_biot_savart_u(1.,j*dy,nx,ny,dx,dy,vort);
   }
   for (j=1;j<ny;j++) {
     psi[nx-1][j] = psi[nx-1][j-1] + 0.5*(vel[j]+vel[j-1])*dy;
     //fprintf(stderr,"y=%g  u=%g  psi=%g\n",j*dy,vel[j],psi[nx-1][j]);
   }

   // finally, march across top row (j=ny-1, y=yf)
   #pragma omp parallel for private(i,thisvel)
   for (i=0;i<nx;i++) {
     find_biot_savart (i*dx,yf,nx,ny,dx,dy,vort,thisvel);
     vel[i] = freestream[1] + thisvel[1];
     //vel[i] = freestream[1] + find_biot_savart_v(i*dx,yf,nx,ny,dx,dy,vort);
   }
   for (i=1;i<nx;i++) {
     psi[i][ny-1] = psi[i-1][ny-1] - 0.5*(vel[i]+vel[i-1])*dx;
     //fprintf(stderr,"x=%g  v=%g  psi=%g\n",i*dx,vel[i],psi[i][ny-1]);
   }
   //exit(0);
   free_1d_array_f(vel);

   return(0);
}

float find_biot_savart_u (float xp,float yp,int nx,int ny,
                        float dx,float dy,float **vort) {
  int i,j;
  float xdist,ydist,distsq;
  float uvel = 0.;

  for (i=0;i<nx;i++) {
    xdist = i*dx - xp;
    for (j=0;j<ny;j++) {
      ydist = j*dy - yp;
      distsq = xdist*xdist+ydist*ydist;
      if (distsq > 1.e-20) {
        uvel += ydist * vort[i][j] / distsq;
      }
    }
  }
  uvel *= dx*dy/(2.*M_PI);

  return uvel;
}

float find_biot_savart_v (float xp,float yp,int nx,int ny,
                        float dx,float dy,float **vort) {
  int i,j;
  float xdist,ydist,distsq;
  float vvel = 0.;

  for (i=0;i<nx;i++) {
    xdist = i*dx - xp;
    for (j=0;j<ny;j++) {
      ydist = j*dy - yp;
      distsq = xdist*xdist+ydist*ydist;
      if (distsq > 1.e-20) {
        vvel -= xdist * vort[i][j] / distsq;
      }
    }
  }
  vvel *= dx*dy/(2.*M_PI);

  return vvel;
}

int find_biot_savart (float xp,float yp,int nx,int ny,
                      float dx,float dy,float **vort,
                      float *vel) {
  int i,j;
  static int istart = 100000;
  static int iend = -1;
  static int jstart = 100000;
  static int jend = -1;
  float xdist,ydist,distsq;

  // initialize some arrays
  if (dx < 0. || dy < 0.) {

    istart = 100000;
    iend = -1;
    jstart = 100000;
    jend = -1;
    
    // which columns have non-zero circulation?
    // which rows have non-zero circulation?
    for (i=0;i<nx;i++) {
      for (j=0;j<ny;j++) {
        if (abs(vort[i][j]) > 1.e-20) {
          if (i<istart) istart=i;
          if (i+1>iend) iend=i+1;
          if (j<jstart) jstart=j;
          if (j+1>jend) jend=j+1;
        }
      }
    }

    //fprintf(stderr,"%d %d  %d %d\n",istart,iend,jstart,jend);
  }

  // now, only iterate over the necessary rows/columns
  vel[0] = 0.;
  vel[1] = 0.;
  for (i=istart;i<iend;i++) {
    xdist = i*dx - xp;
    for (j=jstart;j<jend;j++) {
      ydist = j*dy - yp;
      distsq = xdist*xdist+ydist*ydist;
      if (distsq > 1.e-20) {
        vel[0] += ydist * vort[i][j] / distsq;
        vel[1] -= xdist * vort[i][j] / distsq;
      }
    }
  }
  vel[0] *= dx*dy/(2.*M_PI);
  vel[1] *= dx*dy/(2.*M_PI);

  return ((iend-istart-1)*(jend-jstart-1));
}


/*
 * make a vector field solenoidal by subtracting the irrotational component
 */
int make_solenoidal_3d(int nx,int ny,int nz,int xbdry,int ybdry,int zbdry,
   float ***u,float ***v,float ***w) {

   int maxIter = 100;
   float maxAllowedDiv = 2.;
   int iter,i,j,k;
   int nxm1 = nx-1;
   int nym1 = ny-1;
   int nzm1 = nz-1;
   int nxp = nx;
   int nyp = ny;
   int nzp = nz;
   int bcux,bcvx,bcwx,bcuy,bcvy,bcwy,bcuz,bcvz,bcwz;
   float elmbda = 0.0;
   // these are for both
   int ierr;
   float pert;
   float xs = 0.0;
   float ys = 0.0;
   float zs = 0.0;
   float xf = 1.0;
   float yf = (float)((float)(ny-1)/(float)(nx-1));
   float zf = (float)((float)(nz-1)/(float)(nx-1));
   static int iworksize = 4000000;
   static float *work;	// increasing this does not fix hwscrt's problem
   static int allocated_arrays = FALSE;
   static int must_initialize = TRUE;
   static float ***rhs;
   //static float ***phi;
   static float ***gradphi[3];
   static float **bduxs,**bduxf,**bduys,**bduyf,**bduzs,**bduzf;
   static float **bdvxs,**bdvxf,**bdvys,**bdvyf,**bdvzs,**bdvzf;
   static float **bdwxs,**bdwxf,**bdwys,**bdwyf,**bdwzs,**bdwzf;
   static float **xswx,**xswy,**xswz,**xfwx,**xfwy,**xfwz;
   static float **yswx,**yswy,**yswz,**yfwx,**yfwy,**yfwz;
   static float **zswx,**zswy,**zswz,**zfwx,**zfwy,**zfwz;
   float div,maxdiv;
   char outfileroot[70];

   // fprintf(stderr,"in find_vels_3d\n"); fflush(stderr);

   // allocate memory for the streamfunction, and the boundary
   //    vorticity arrays
   if (!allocated_arrays) {

      // always allocate work array
      work = allocate_1d_array_f(iworksize);

      // allocate rhs
      rhs = allocate_3d_array_f(nx,ny,nz);

      // allocate for unknown phi
      //phi = allocate_3d_array_f(nx,ny,nz);
      gradphi[0] = allocate_3d_array_f(nx,ny,nz);
      gradphi[1] = allocate_3d_array_f(nx,ny,nz);
      gradphi[2] = allocate_3d_array_f(nx,ny,nz);

      // allocate memory for the boundary condition arrays
      // if (xbdry == WALL) {
         bduxs = allocate_2d_array_f(ny,nz);
         bduxf = allocate_2d_array_f(ny,nz);
      // }
      // if (ybdry == WALL) {
         bduys = allocate_2d_array_f(nx,nz);
         bduyf = allocate_2d_array_f(nx,nz);
      // }
      // if (zbdry == WALL) {
         bduzs = allocate_2d_array_f(nx,ny);
         bduzf = allocate_2d_array_f(nx,ny);
      // }
      bdvxs = allocate_2d_array_f(ny,nz);
      bdvxf = allocate_2d_array_f(ny,nz);
      bdvys = allocate_2d_array_f(nx,nz);
      bdvyf = allocate_2d_array_f(nx,nz);
      bdvzs = allocate_2d_array_f(nx,ny);
      bdvzf = allocate_2d_array_f(nx,ny);
      bdwxs = allocate_2d_array_f(ny,nz);
      bdwxf = allocate_2d_array_f(ny,nz);
      bdwys = allocate_2d_array_f(nx,nz);
      bdwyf = allocate_2d_array_f(nx,nz);
      bdwzs = allocate_2d_array_f(nx,ny);
      bdwzf = allocate_2d_array_f(nx,ny);

      xswx = allocate_2d_array_f(ny,nz);
      xswy = allocate_2d_array_f(ny,nz);
      xswz = allocate_2d_array_f(ny,nz);
      xfwx = allocate_2d_array_f(ny,nz);
      xfwy = allocate_2d_array_f(ny,nz);
      xfwz = allocate_2d_array_f(ny,nz);

      yswx = allocate_2d_array_f(nx,nz);
      yswy = allocate_2d_array_f(nx,nz);
      yswz = allocate_2d_array_f(nx,nz);
      yfwx = allocate_2d_array_f(nx,nz);
      yfwy = allocate_2d_array_f(nx,nz);
      yfwz = allocate_2d_array_f(nx,nz);

      zswx = allocate_2d_array_f(nx,ny);
      zswy = allocate_2d_array_f(nx,ny);
      zswz = allocate_2d_array_f(nx,ny);
      zfwx = allocate_2d_array_f(nx,ny);
      zfwy = allocate_2d_array_f(nx,ny);
      zfwz = allocate_2d_array_f(nx,ny);

      allocated_arrays = TRUE;
   }


   // iterate this whole thing a few times

   maxdiv = 2.*maxAllowedDiv;
   iter = 0;
   while (iter < maxIter && maxdiv > maxAllowedDiv) {

   iter++;

   // initialize the right-hand-side as the divergence of the vector field
   maxdiv = compute_divergence_3d(nx,ny,nz,xbdry,ybdry,zbdry,u,v,w,rhs);

   // Use fishpak FFT solver ----------------------------------------------

   // for all wall boundaries, set zero thru-flow and derivative
   //if (slip_wall) {
   if (TRUE) {
      // allow slip at walls (LES-like)
      if (xbdry == WALL) {
         for (j=0;j<ny;j++) {
            for (k=0;k<nz;k++) {
               u[0][j][k] = 0.0;
               u[nxm1][j][k] = 0.0;
               bdvxs[j][k] = 0.0;
               bdvxf[j][k] = 0.0;
               bdwxs[j][k] = 0.0;
               bdwxf[j][k] = 0.0;
            }
         }
      }
      if (ybdry == WALL) {
         for (i=0;i<nx;i++) {
            for (k=0;k<nz;k++) {
               v[i][0][k] = 0.0;
               v[i][nym1][k] = 0.0;
               bduys[i][k] = 0.0;
               bduyf[i][k] = 0.0;
               bdwys[i][k] = 0.0;
               bdwyf[i][k] = 0.0;
            }
         }
      }
      if (zbdry == WALL) {
         for (i=0;i<nx;i++) {
            for (j=0;j<ny;j++) {
               w[i][j][0] = 0.0;
               w[i][j][nzm1] = 0.0;
               bduzs[i][j] = 0.0;
               bduzf[i][j] = 0.0;
               bdvzs[i][j] = 0.0;
               bdvzf[i][j] = 0.0;
            }
         }
      }
   } else {
      // no-slip wall (viscous simulation)
      if (xbdry == WALL) {
         for (j=0;j<ny;j++) {
            for (k=0;k<nz;k++) {
               u[0][j][k] = 0.0;
               u[nxm1][j][k] = 0.0;
               // and tangential components should be zero, too
               v[0][j][k] = 0.0;
               v[nxm1][j][k] = 0.0;
               w[0][j][k] = 0.0;
               w[nxm1][j][k] = 0.0;
            }
         }
      }
      if (ybdry == WALL) {
         for (i=0;i<nx;i++) {
            for (k=0;k<nz;k++) {
               v[i][0][k] = 0.0;
               v[i][nym1][k] = 0.0;
               // and tangential components should be zero, too
               u[i][0][k] = 0.0;
               u[i][nym1][k] = 0.0;
               w[i][0][k] = 0.0;
               w[i][nym1][k] = 0.0;
            }
         }
      }
      if (zbdry == WALL) {
         for (i=0;i<nx;i++) {
            for (j=0;j<ny;j++) {
               w[i][j][0] = 0.0;
               w[i][j][nzm1] = 0.0;
               // and tangential components should be zero, too
               u[i][j][0] = 0.0;
               u[i][j][nzm1] = 0.0;
               v[i][j][0] = 0.0;
               v[i][j][nzm1] = 0.0;
            }
         }
      }
   }

   // write out RHS
   // sprintf(outfileroot,"dcvx_%04d",step);
   // write_output_3d(outfileroot,nx,ny,nz,u,-500,1000,2,FALSE);
   // sprintf(outfileroot,"dcvy_%04d",step);
   // write_output_3d(outfileroot,nx,ny,nz,v,-500,1000,2,FALSE);
   // sprintf(outfileroot,"dcvz_%04d",step);
   // write_output_3d(outfileroot,nx,ny,nz,w,-500,1000,2,FALSE);

   // set bc[uvw][xyz] flag boundary conditions
   // for '1' (value specified at boundary), the value is in u/v/w,
   //   and bd[xyz][sf] is not even used (no-slip walls)
   // for '3' (derivative specified at bdry), u/v/w is the rhs, and
   //   bd[xyz][sf] is the derivative (slip wall case)
   if (xbdry == PERIODIC) {
      bcux = 0;
      bcvx = 0;
      bcwx = 0;
   } else if (xbdry == WALL) {
      /* both sides are walls */
      //if (slip_wall) {
      if (TRUE) {
         bcux = 1;
         bcvx = 3;
         bcwx = 3;
      } else {
         bcux = 1;
         bcvx = 1;
         bcwx = 1;
      }
   } else {
      fprintf(stdout,"ERROR (find_vels_3d): BCs for x-dir invalid\n");
      fprintf(stdout,"exiting.\n");
      exit(0);
   }
   if (ybdry == PERIODIC) {
      bcuy = 0;
      bcvy = 0;
      bcwy = 0;
   } else if (ybdry == WALL) {
      /* both sides are walls */
      //if (slip_wall) {
      if (TRUE) {
         bcuy = 3;
         bcvy = 1;
         bcwy = 3;
      } else {
         bcuy = 1;
         bcvy = 1;
         bcwy = 1;
      }
   } else {
      fprintf(stdout,"ERROR (find_vels_3d): BCs for y-dir invalid\n");
      fprintf(stdout,"exiting.\n");
      exit(0);
   }
   if (zbdry == PERIODIC) {
      bcuz = 0;
      bcvz = 0;
      bcwz = 0;
   } else if (zbdry == WALL) {
      /* both sides are walls */
      //if (slip_wall) {
      if (TRUE) {
         bcuz = 3;
         bcvz = 3;
         bcwz = 1;
      } else {
         bcuz = 1;
         bcvz = 1;
         bcwz = 1;
      }
   } else {
      fprintf(stdout,"ERROR (find_vels_3d): BCs for z-dir invalid\n");
      fprintf(stdout,"exiting.\n");
      exit(0);
   }

   // make the call to the external solver, note that x and y are backwards...

#ifdef IS_3D
   hw3crt_(&zs,&zf,&nzm1,&bcuz,bduxs[0],bduxf[0],
           &ys,&yf,&nym1,&bcuy,bduys[0],bduyf[0],
           &xs,&xf,&nxm1,&bcux,bduzs[0],bduzf[0],
           &elmbda,&nzp,&nyp,rhs[0][0],&pert,&ierr,work);
#endif

   //fprintf(stdout,"  hwscrt used %d locations in work array\n",(int)(work[0]));
   //fprintf(stdout,"    err=%d, pert=%g\n",ierr,pert);

   // find gradients of solution (phi)
   find_gradient_of_scalar_3d(nx,ny,nz,xbdry,ybdry,zbdry,rhs,gradphi);
   //sprintf(outfileroot,"gradphi_%04d",0);
   //write_3d_vtk(outfileroot,nx,ny,nz,gradphi[0],gradphi[1],gradphi[2]);

   // subtract them off of the input vector field
   for (i=0; i<nx; i++)
   for (j=0; j<ny; j++)
   for (k=0; k<nz; k++) {
     //u[i][j][k] -= 1.0*gradphi[0][i][j][k];
     //v[i][j][k] -= 1.0*gradphi[1][i][j][k];
     //w[i][j][k] -= 1.0*gradphi[2][i][j][k];
     // over-relax!
     u[i][j][k] -= 1.5*gradphi[0][i][j][k];
     v[i][j][k] -= 1.5*gradphi[1][i][j][k];
     w[i][j][k] -= 1.5*gradphi[2][i][j][k];
   }
   //sprintf(outfileroot,"newvort_%04d",0);
   //write_3d_vtk(outfileroot,nx,ny,nz,u,v,w);

   }  // end loop over iter
   printf("\n");
   fflush(stdout);

   // now, out of curiosity, what is the divergence?
   //compute_divergence_3d(nx,ny,nz,xbdry,ybdry,zbdry,u,v,w,rhs);
   //exit(0);

   return(0);
}


/*
 * differentiate the a scalar to find the velocity
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
 */
int find_gradient_of_scalar_2nd_2d (int nx,int ny,int xbdry,int ybdry,float **psi,float **u,float umult,float **v,float vmult) {

   int i,j;
   int nxm1 = nx-1;
   int nym1 = ny-1;
   FLOAT hxi,hyi;
   hxi = 0.5*umult*(nxm1)/SCALE;	// negative because vel is neg gradient
   hyi = 0.5*vmult*(nxm1);		// and *nxm1 because xsize=1

   // do the middle first
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

   return(0);
}


/*
 * differentiate the a scalar to find the velocity - 3D now
 *
 * done: use a new method that is 2nd order everywhere
 *
 * psi is input
 * u is output
 */
int find_gradient_of_scalar_3d (int nx,int ny,int nz,int xbdry,int ybdry,int zbdry,
   float ***psi,float ****u) {

   int i,j,k;
   int nxm1 = nx-1;
   int nym1 = ny-1;
   int nzm1 = nz-1;
   FLOAT hxi,hyi,hzi;
   hxi = 0.5*(nxm1);		// negative because vel is neg gradient
   hyi = 0.5*(nxm1);		// and *nxm1 because xsize=1
   hzi = 0.5*(nxm1);

   // partial in x-direction
   // do the middle first
   for (i=1; i<nxm1; i++)
   for (j=0; j<ny; j++)
   for (k=0; k<nz; k++) {
      u[0][i][j][k] = (psi[i+1][j][k]-psi[i-1][j][k])*hxi;
   }
   // now, do left and right sides (-x and +x)
   if (xbdry == WALL) {
      i = 0;
      for (j=0; j<ny; j++)
      for (k=0; k<nz; k++) {
         u[0][i][j][k] = (-1.*psi[i+2][j][k]+4.*psi[i+1][j][k]-3.*psi[i][j][k])*hxi;
      }
      i = nxm1;
      for (j=0; j<ny; j++)
      for (k=0; k<nz; k++) {
         u[0][i][j][k] = (psi[i-2][j][k]-4.*psi[i-1][j][k]+3.*psi[i][j][k])*hxi;
      }
   } else if (xbdry == PERIODIC) {
      i = 0;
      for (j=0; j<ny; j++)
      for (k=0; k<nz; k++) {
         u[0][i][j][k] = (psi[i+1][j][k]-psi[nxm1-1][j][k])*hxi;
         u[0][nxm1][j][k] = u[0][i][j][k];
      }
   }
 
   // partial in y-direction
   for (i=0; i<nx; i++)
   for (j=1; j<nym1; j++)
   for (k=0; k<nz; k++) {
      u[1][i][j][k] = (psi[i][j+1][k]-psi[i][j-1][k])*hyi;
   }
   // now, do front and back sides (-y and +y)
   if (ybdry == WALL) {
      j = 0;
      for (i=0; i<nx; i++)
      for (k=0; k<nz; k++) {
         u[1][i][j][k] = (-1.*psi[i][j+2][k]+4.*psi[i][j+1][k]-3.*psi[i][j][k])*hyi;
      }
      j = nym1;
      for (i=0; i<nx; i++)
      for (k=0; k<nz; k++) {
         u[1][i][j][k] = (psi[i][j-2][k]-4.*psi[i][j-1][k]+3.*psi[i][j][k])*hyi;
      }
   } else if (ybdry == PERIODIC) {
      j = 0;
      for (i=0; i<nx; i++)
      for (k=0; k<nz; k++) {
         u[1][i][j][k] = (psi[i][j+1][k]-psi[i][nym1-1][k])*hyi;
         u[1][i][nym1][k] = u[1][i][j][k];
      }
   }

   // partial in z-direction
   for (i=0; i<nx; i++)
   for (j=0; j<ny; j++)
   for (k=1; k<nzm1; k++) {
      u[2][i][j][k] = (psi[i][j][k+1]-psi[i][j][k-1])*hzi;
   }
   // now, do bottom and top sides (-z and +z)
   if (zbdry == WALL) {
      k = 0;
      for (i=0; i<nx; i++)
      for (j=0; j<ny; j++) {
         u[2][i][j][k] = (-1.*psi[i][j][k+2]+4.*psi[i][j][k+1]-3.*psi[i][j][k])*hzi;
      }
      k = nzm1;
      for (i=0; i<nx; i++)
      for (j=0; j<ny; j++) {
         u[2][i][j][k] = (psi[i][j][k-2]-4.*psi[i][j][k-1]+3.*psi[i][j][k])*hzi;
      }
   } else if (zbdry == PERIODIC) {
      k = 0;
      for (i=0; i<nx; i++)
      for (j=0; j<ny; j++) {
         u[2][i][j][k] = (psi[i][j][k+1]-psi[i][j][nzm1-1])*hzi;
         u[2][i][j][nzm1] = u[2][i][j][k];
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
   hxi = 0.5*(nxm1)/SCALE;
   hyi = 0.5*(nxm1);		// and *nxm1 because xsize=1

   // do the middle first
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
 * moc_advect_2d is the implicit scheme for updating the vorticity values
 */
int moc_advect_2d (int nx,int ny,int xbdry,int ybdry,float **u,float **v,
   float ***in,float ***out,int sc_cnt,float dt) {

   //int useCIC = FALSE;
   //int useCIC = TRUE;
   enum interpMeth {
      cic, tsc, m4p } interp = m4p;
   int useOneStep = FALSE;
   //int useOneStep = TRUE;
   int i,j,k;
   float yf = (float)(ny-1)/(float)(nx-1);
   float px,py,newx,newy,temp;
   float u0,u1,v0,v1;
   float outvals[MAX_SCALARS];

   // fprintf(stderr,"in moc_advect_2d\n"); fflush(stderr);

   // for each point, march backwards in space and time to find the value
   if (useOneStep) {
      // use a 1-step method
      for (i=0;i<nx;i++) {
         px = (float)i/(float)(nx-1);
         for (j=0;j<ny;j++) {
            // py = (float)j/(float)(ny-1);
            py = (float)j/(float)(nx-1);
            newx = px-dt*u[i][j];
            newy = py-dt*v[i][j];
            // 2007-01-20 comment these out
            //if (newx > 1.0-EPSILON) newx = 1.0-EPSILON;
            //if (newx < 0.0+EPSILON) newx = EPSILON;
            //if (newy > yf-EPSILON) newy = yf-EPSILON;
            //if (newy < 0.0+EPSILON) newy = EPSILON;
            // temp = interpolate_using_area_2d(nx,ny,xbdry,ybdry,in,newx,newy);
            // temp = interpolate_using_M4p_2d(nx,ny,xbdry,ybdry,in,newx,newy);
            if (interp == cic)
               interpolate_array_using_CIC_2d(nx,ny,xbdry,ybdry,in,newx,newy,sc_cnt,outvals);
            else if (interp == tsc)
               interpolate_array_using_TSC_2d(nx,ny,xbdry,ybdry,in,newx,newy,sc_cnt,outvals);
            else
               interpolate_array_using_M4p_2d(nx,ny,xbdry,ybdry,in,newx,newy,sc_cnt,outvals);
            // out[i][j] = temp;
            for (k=0; k<sc_cnt; k++) out[k][i][j] = outvals[k];
         }
      }
   } else {
      // use a two-step method
      #pragma omp parallel for private(i,px,j,py,u0,u1,v0,v1,newx,newy,k,outvals)
      for (i=0;i<nx;i++) {
         px = (float)i/(float)(nx-1);
         for (j=0;j<ny;j++) {
            // py = (float)j/(float)(ny-1);
            py = (float)j/(float)(nx-1);
            // find vel at px,py
            u0 = u[i][j];
            v0 = v[i][j];
            // find position 1 explicit Euler step backwards
            newx = px-dt*u0;
            if (xbdry == WALL || xbdry == OPEN) {
               if (newx > 1.0-EPSILON) newx = 1.0-EPSILON;
               if (newx < 0.0+EPSILON) newx = EPSILON;
            }
            newy = py-dt*v0;
            if (ybdry == WALL || ybdry == OPEN) {
               if (newy > yf-EPSILON) newy = yf-EPSILON;
               if (newy < 0.0+EPSILON) newy = EPSILON;
            }
            // find vel at newx,newy
            if (interp == cic)
               (void) interpolate_vel_using_CIC_2d(nx,ny,xbdry,ybdry,u,v,newx,newy,&u1,&v1);
            else if (interp == tsc)
               (void) interpolate_vel_using_TSC_2d(nx,ny,xbdry,ybdry,u,v,newx,newy,&u1,&v1);
            else
               (void) interpolate_vel_using_M4p_2d(nx,ny,xbdry,ybdry,u,v,newx,newy,&u1,&v1);
            // find position back 1 step using average of two velocities
            newx = px-dt*0.5*(u0+u1);
            if (xbdry == WALL || xbdry == OPEN) {
               if (newx > 1.0-EPSILON) newx = 1.0-EPSILON;
               if (newx < 0.0+EPSILON) newx = EPSILON;
            }
            newy = py-dt*0.5*(v0+v1);
            if (ybdry == WALL || ybdry == OPEN) {
               if (newy > yf-EPSILON) newy = yf-EPSILON;
               if (newy < 0.0+EPSILON) newy = EPSILON;
            }
            // finally, return scalar-valued vorticity from there
            if (interp == cic)
               interpolate_array_using_CIC_2d(nx,ny,xbdry,ybdry,in,newx,newy,sc_cnt,outvals);
            else if (interp == tsc)
               interpolate_array_using_TSC_2d(nx,ny,xbdry,ybdry,in,newx,newy,sc_cnt,outvals);
            else
               interpolate_array_using_M4p_2d(nx,ny,xbdry,ybdry,in,newx,newy,sc_cnt,outvals);
            // out[i][j] = temp;
            for (k=0; k<sc_cnt; k++) out[k][i][j] = outvals[k];
         }
      }
   }

   return(0);
}


/*
 * moc_advect_3d is the implicit scheme for updating the vorticity values
 */
int moc_advect_3d (int nx,int ny,int nz,int xbdry,int ybdry,int zbdry,
   float ****u,float ****in,float ****out,int sc_cnt,float dt) {

   int i,j,k,l;
   float yf = (float)(ny-1)/(float)(nx-1);
   float zf = (float)(nz-1)/(float)(nx-1);
   float px,py,pz,newx,newy,newz,temp;
   float u0,u1,v0,v1,w0,w1;
   float outvels[MAX_SCALARS];
   float outvals[MAX_SCALARS];

   // fprintf(stderr,"in moc_advect_3d\n"); fflush(stderr);

   // for each point, march backwards in space and time to find the value
   // use a two-step method
   for (i=0;i<nx;i++) {
      px = (float)i/(float)(nx-1);
      for (j=0;j<ny;j++) {
         py = (float)j/(float)(nx-1);
         for (k=0;k<nz;k++) {
            pz = (float)k/(float)(nx-1);

            // find vel at px,py,pz
            u0 = u[XV][i][j][k];
            v0 = u[YV][i][j][k];
            w0 = u[ZV][i][j][k];

            // find position 1 explicit Euler step backwards
            newx = px-dt*u0;
            if (xbdry == WALL) {
               if (newx > 1.0-EPSILON) newx = 1.0-EPSILON;
               if (newx < 0.0+EPSILON) newx = EPSILON;
            }
            newy = py-dt*v0;
            if (ybdry == WALL) {
               if (newy > yf-EPSILON) newy = yf-EPSILON;
               if (newy < 0.0+EPSILON) newy = EPSILON;
            }
            newz = pz-dt*w0;
            if (zbdry == WALL) {
               if (newz > zf-EPSILON) newz = zf-EPSILON;
               if (newz < 0.0+EPSILON) newz = EPSILON;
            }

            // find vel at newx,newy
            //interpolate_array_using_CIC_3d(nx,ny,nz,xbdry,ybdry,zbdry,u,newx,newy,newz,3,outvels);
            interpolate_array_using_M4p_3d(nx,ny,nz,xbdry,ybdry,zbdry,u,newx,newy,newz,3,outvels);
            u1 = outvels[0];
            v1 = outvels[1];
            w1 = outvels[2];

            // find position back 1 step using average of two velocities
            newx = px-dt*0.5*(u0+u1);
            if (xbdry == WALL) {
               if (newx > 1.0-EPSILON) newx = 1.0-EPSILON;
               if (newx < 0.0+EPSILON) newx = EPSILON;
            }
            newy = py-dt*0.5*(v0+v1);
            if (ybdry == WALL) {
               if (newy > yf-EPSILON) newy = yf-EPSILON;
               if (newy < 0.0+EPSILON) newy = EPSILON;
            }
            newz = pz-dt*0.5*(w0+w1);
            if (zbdry == WALL) {
               if (newz > zf-EPSILON) newz = zf-EPSILON;
               if (newz < 0.0+EPSILON) newz = EPSILON;
            }

            // finally, return scalar-valued vorticity from there
            //interpolate_array_using_CIC_3d(nx,ny,nz,xbdry,ybdry,zbdry,in,newx,newy,newz,sc_cnt,outvals);
            interpolate_array_using_M4p_3d(nx,ny,nz,xbdry,ybdry,zbdry,in,newx,newy,newz,sc_cnt,outvals);

            // add those to the output
            for (l=0; l<sc_cnt; l++) out[l][i][j][k] = outvals[l];
         }
      }
   }

   return(0);
}


/*
 * advect the particles using an explicit scheme
 *
 * for this first step, just assume the particles move
 * with the flow---no ballistic coefficient
 */
int explicit_particle_move_3d(int nx,int ny,int nz,int xbdry,int ybdry,int zbdry,
   float ****u,float dt,float re,int pnum,float **ploc,float **pvel) {

   int i;
   float temploc[3],tempvel[3];
   float dr,len,rvec[3];

   // first, apply some random motion, scaled by sqrt(dt) and 1/Re
   if (TRUE) {
      for (i=0; i<pnum; i++) {
         // find a random direction
         len = 2.;
         while (len > 1.) {
            rvec[0] = 2.*rand()/(RAND_MAX+1.0) - 1.;
            rvec[1] = 2.*rand()/(RAND_MAX+1.0) - 1.;
            rvec[2] = 2.*rand()/(RAND_MAX+1.0) - 1.;
            len = sqrt(rvec[0]*rvec[0]+rvec[1]*rvec[1]+rvec[2]*rvec[2]);
         }
         // scale it by a proper distance (must use Gaussian!)
         // dr = pow(rand()/(RAND_MAX+1.0),2)*sqrt(dt/re);
         dr = (float)(gaussian(0.0,(double)(sqrt(dt/re))));
         // add it to the current locations
         ploc[i][0] += dr*rvec[0]/len;
         ploc[i][1] += dr*rvec[1]/len;
         ploc[i][2] += dr*rvec[2]/len;
         // and make sure it didn't go out of bounds
         if (xbdry == WALL) {
            if (ploc[i][0] < 0.0) ploc[i][0] = -ploc[i][0];
            if (ploc[i][0] > 1.0) ploc[i][0] = 2.0-ploc[i][0];
         }
         if (ybdry == WALL) {
            if (ploc[i][1] < 0.0) ploc[i][1] = -ploc[i][1];
            if (ploc[i][1] > 1.0) ploc[i][1] = 2.0-ploc[i][1];
         }
         if (zbdry == WALL) {
            if (ploc[i][2] < 0.0) ploc[i][2] = -ploc[i][2];
            if (ploc[i][2] > 1.0) ploc[i][2] = 2.0-ploc[i][2];
         }
      }
   }

   if (FALSE) {
      // find new particle velocities and locations, Euler-style
      for (i=0; i<pnum; i++) {
         //interpolate_array_using_M4p_3d(nx,ny,nz,xbdry,ybdry,zbdry,u,
         interpolate_array_using_CIC_3d(nx,ny,nz,xbdry,ybdry,zbdry,u,
            ploc[i][0],ploc[i][1],ploc[i][2],3,pvel[i]);
         ploc[i][0] += dt*pvel[i][0];
         ploc[i][1] += dt*pvel[i][1];
         ploc[i][2] += dt*pvel[i][2];
      }
   } else {
      // use 2-step Runge-Kutta
      for (i=0; i<pnum; i++) {
         //interpolate_array_using_M4p_3d(nx,ny,nz,xbdry,ybdry,zbdry,u,
         interpolate_array_using_CIC_3d(nx,ny,nz,xbdry,ybdry,zbdry,u,
            ploc[i][0],ploc[i][1],ploc[i][2],3,tempvel);
         temploc[0] = ploc[i][0] + dt*tempvel[0];
         temploc[1] = ploc[i][1] + dt*tempvel[1];
         temploc[2] = ploc[i][2] + dt*tempvel[2];
         //interpolate_array_using_M4p_3d(nx,ny,nz,xbdry,ybdry,zbdry,u,
         interpolate_array_using_CIC_3d(nx,ny,nz,xbdry,ybdry,zbdry,u,
            temploc[0],temploc[1],temploc[2],3,pvel[i]);
         pvel[i][0] = 0.5*(pvel[i][0]+tempvel[0]);
         pvel[i][1] = 0.5*(pvel[i][1]+tempvel[1]);
         pvel[i][2] = 0.5*(pvel[i][2]+tempvel[2]);
         ploc[i][0] += dt*pvel[i][0];
         ploc[i][1] += dt*pvel[i][1];
         ploc[i][2] += dt*pvel[i][2];
      }
   }

   return(0);
}


/* 
 * compute the divergence of a vector field, globally 2nd order
 */
float compute_divergence_3d(int nx,int ny,int nz,int xbdry,int ybdry,int zbdry,
   float ***u,float ***v,float ***w,float ***out) {

   int i,j,k;
   int nxm1 = nx-1;
   int nym1 = ny-1;
   int nzm1 = nz-1;
   int xs,xf,ys,yf,zs,zf;
   float mult = 0.5*nxm1;
   float maxdiv = 0.;

   // do the middle part of the field, regardless of the BCs
   for (i=1;i<nxm1;i++) {
     for (j=1;j<nym1;j++) {
       for (k=1;k<nzm1;k++) {
         out[i][j][k] = mult*(
           +u[i+1][j][k]-u[i-1][j][k]
           +v[i][j+1][k]-v[i][j-1][k]
           +w[i][j][k+1]-w[i][j][k-1]);
       }
     }
   }

   // now, do the boundaries, faces first (one bdry variable)
   if (xbdry == WALL) {
     fprintf(stderr,"ERROR (compute_divergence_3d): not programmed\n");
     exit(0);
     for (j=1;j<nym1;j++) {
     for (k=1;k<nzm1;k++) {
       //out[0][j][k] = in[0][j][k] + mult*(
       //out[nxm1][j][k] = in[nxm1][j][k] + mult*(
     }
     }
   } else if (xbdry == PERIODIC) {
     for (j=1;j<nym1;j++) {
       for (k=1;k<nzm1;k++) {
         out[0][j][k] = mult*(
           +u[1][j][k]-u[nxm1-1][j][k]
           +v[0][j+1][k]-v[0][j-1][k]
           +w[0][j][k+1]-w[0][j][k-1]);
         out[nxm1][j][k] = out[0][j][k];
       }
     }
   }
   if (ybdry == WALL) {
     fprintf(stderr,"ERROR (compute_divergence_3d): not programmed\n");
     exit(0);
     for (i=1;i<nxm1;i++) {
     for (k=1;k<nzm1;k++) {
       //out[i][0][k] = in[i][0][k] + mult*(
       //out[i][nym1][k] = in[i][nym1][k] + mult*(
     }
     }
   } else if (ybdry == PERIODIC) {
     for (i=1;i<nxm1;i++) {
       for (k=1;k<nzm1;k++) {
         out[i][0][k] = mult*(
           +u[i+1][0][k]-u[i-1][0][k]
           +v[i][1][k]-v[i][nym1-1][k]
           +w[i][0][k+1]-w[i][0][k-1]);
         out[i][nym1][k] = out[i][0][k];
       }
     }
   }
   if (zbdry == WALL) {
     fprintf(stderr,"ERROR (compute_divergence_3d): not programmed\n");
     exit(0);
     for (i=1;i<nxm1;i++) {
     for (j=1;j<nym1;j++) {
       //out[i][j][0] = in[i][j][0] + mult*(
       //out[i][j][nzm1] = in[i][j][nzm1] + mult*(
     }
     }
   } else if (zbdry == PERIODIC) {
     for (i=1;i<nxm1;i++) {
       for (j=1;j<nym1;j++) {
         out[i][j][0] = mult*(
           +u[i+1][j][0]-u[i-1][j][0]
           +v[i][j+1][0]-v[i][j-1][0]
           +w[i][j][1]-w[i][j][nzm1-1]);
         out[i][j][nzm1] = out[i][j][0];
       }
     }
   }

   // now, the 12 edges! ouch; first the edges in which x varies
   if (ybdry == WALL && zbdry == WALL) {
     fprintf(stderr,"ERROR (compute_divergence_3d): not programmed\n");
     exit(0);
     for (i=1;i<nxm1;i++) {
       //out[i][0][0] = in[i][0][0] + mult*(
       //out[i][0][nzm1] = in[i][0][nzm1] + mult*(
       //out[i][nym1][nzm1] = in[i][nym1][nzm1] + mult*(
       //out[i][nym1][0] = in[i][nym1][0] + mult*(
     }
   } else if (ybdry == PERIODIC && zbdry == WALL) {
     fprintf(stderr,"ERROR (compute_divergence_3d): not programmed\n");
     exit(0);
   } else if (ybdry == WALL && zbdry == PERIODIC) {
     fprintf(stderr,"ERROR (compute_divergence_3d): not programmed\n");
     exit(0);
   } else if (ybdry == PERIODIC && zbdry == PERIODIC) {
     //fprintf(stderr,"ERROR (compute_divergence_3d): not programmed\n");
     //exit(0);
     for (i=1;i<nxm1;i++) {
       out[i][0][0] = mult*(
         +u[i+1][0][0]-u[i-1][0][0]
         +v[i][1][0]-v[i][nym1-1][0]
         +w[i][0][1]-w[i][0][nzm1-1]);
       out[i][0][nzm1] = out[i][0][0];
       out[i][nym1][nzm1] = out[i][0][0];
       out[i][nym1][0] = out[i][0][0];
     }
   }
   // then the edges for which y varies
   if (xbdry == WALL && zbdry == WALL) {
     fprintf(stderr,"ERROR (compute_divergence_3d): not programmed\n");
     exit(0);
     for (j=1;j<nym1;j++) {
       //out[0][j][0] = in[0][j][0] + mult*(
       //out[0][j][nzm1] = in[0][j][nzm1] + mult*(
       //out[nxm1][j][nzm1] = in[nxm1][j][nzm1] + mult*(
       //out[nxm1][j][0] = in[nxm1][j][0] + mult*(
     }
   } else if (ybdry == PERIODIC && zbdry == WALL) {
     fprintf(stderr,"ERROR (compute_divergence_3d): not programmed\n");
     exit(0);
   } else if (ybdry == WALL && zbdry == PERIODIC) {
     fprintf(stderr,"ERROR (compute_divergence_3d): not programmed\n");
     exit(0);
   } else if (ybdry == PERIODIC && zbdry == PERIODIC) {
     //fprintf(stderr,"ERROR (compute_divergence_3d): not programmed\n");
     //exit(0);
     for (j=1;j<nym1;j++) {
       out[0][j][0] = mult*(
         +u[1][j][0]-u[nxm1-1][j][0]
         +v[0][j+1][0]-v[0][j-1][0]
         +w[0][j][1]-w[0][j][nzm1-1]);
       out[0][j][nzm1] = out[0][j][0];
       out[nxm1][j][nzm1] = out[0][j][0];
       out[nxm1][j][0] = out[0][j][0];
     }
   }
   // then the edges for which z varies
   if (xbdry == WALL && ybdry == WALL) {
     fprintf(stderr,"ERROR (compute_divergence_3d): not programmed\n");
     exit(0);
     for (k=1;k<nzm1;k++) {
       //out[0][0][k] = in[0][0][k] + mult*(
       //out[0][nym1][k] = in[0][nym1][k] + mult*(
       //out[nxm1][nym1][k] = in[nxm1][nym1][k] + mult*(
       //out[nxm1][0][k] = in[nxm1][0][k] + mult*(
     }
   } else if (ybdry == PERIODIC && zbdry == WALL) {
     fprintf(stderr,"ERROR (compute_divergence_3d): not programmed\n");
     exit(0);
   } else if (ybdry == WALL && zbdry == PERIODIC) {
     fprintf(stderr,"ERROR (compute_divergence_3d): not programmed\n");
     exit(0);
   } else if (ybdry == PERIODIC && zbdry == PERIODIC) {
     //fprintf(stderr,"ERROR (compute_divergence_3d): not programmed\n");
     //exit(0);
     for (k=1;k<nzm1;k++) {
       out[0][0][k] = mult*(
         +u[1][0][k]-u[nxm1-1][0][k]
         +v[0][1][k]-v[0][nym1-1][k]
         +w[0][0][k+1]-w[0][0][k-1]);
       out[0][nym1][k] = out[0][0][k];
       out[nxm1][nym1][k] = out[0][0][k];
       out[nxm1][0][k] = out[0][0][k];
     }
   }

   // and, finally, the corners
   if (xbdry == PERIODIC && ybdry == PERIODIC && zbdry == PERIODIC) {
      out[0][0][0] = mult*(
         +u[1][0][0]-u[nxm1-1][0][0]
         +v[0][1][0]-v[0][nym1-1][0]
         +w[0][0][1]-w[0][0][nzm1-1]);
      out[0][0][nzm1] = out[0][0][0];
      out[0][nym1][0] = out[0][0][0];
      out[0][nym1][nzm1] = out[0][0][0];
      out[nxm1][0][0] = out[0][0][0];
      out[nxm1][0][nzm1] = out[0][0][0];
      out[nxm1][nym1][0] = out[0][0][0];
      out[nxm1][nym1][nzm1] = out[0][0][0];
   } else {
     fprintf(stderr,"ERROR (compute_divergence_3d): not programmed\n");
     exit(0);
   }

   // find max divergence
   maxdiv = 0;
   for (i=1;i<nxm1;i++) {
     for (j=1;j<nym1;j++) {
       for (k=1;k<nzm1;k++) {
         if (fabs(out[i][j][k]) > maxdiv) {
           maxdiv = fabs(out[i][j][k]);
           //printf("  div %g\n",out[i][j][k]);
           //printf("    at %d %d %d\n",i,j,k);
           //printf("    because of %g %g  %g %g  %g %g\n",u[i+1][j][k],u[i-1][j][k],v[i][j+1][k],v[i][j-1][k],w[i][j][k+1],w[i][j][k-1]);
         }
       }
     }
   }
   //printf("  max div in div3d %g\n",maxdiv);
   printf(".");
   fflush(stdout);
   //exit(0);

   return(maxdiv);
}


/*
 * compute curl of the input, put it in the output arrays
 */
int compute_curl_3d(int nx,int ny,int nz,int xbdry,int ybdry,int zbdry,
   float ***x_in,float ***y_in,float ***z_in,
   float ***x_out,float ***y_out,float ***z_out) {

   int nxm1 = nx-1;
   int nym1 = ny-1;
   int nzm1 = nz-1;
   int i,j,k,im1,ip1,jm1,jp1,km1,kp1;
   int i0im1,i0cx,inip1,incx;
   int j0jm1,j0cy,jnjp1,jncy;
   int k0km1,k0cz,knkp1,kncz;
   float hxi,hyi,hzi,cx,cy,cz;

   // one over cell size, since xsize=1, so this is always nxm1
   hxi = 0.5*(float)(nxm1);
   hyi = 0.5*(float)(nxm1);
   hzi = 0.5*(float)(nxm1);

   /* now, account for periodic BCs by setting some indexes */
   if (xbdry == WALL) {
      i0im1 = 0;
      i0cx = 2.0;
      inip1 = nxm1;
      incx = 2.0;
   } else if (xbdry == PERIODIC) {
      i0im1 = nxm1 - 1;
      i0cx = 1.0;
      inip1 = 1;
      incx = 1.0;
   }

   if (ybdry == WALL) {
      j0jm1 = 0;
      j0cy = 2.0;
      jnjp1 = nym1;
      jncy = 2.0;
   } else if (ybdry == PERIODIC) {
      j0jm1 = nym1 - 1;
      j0cy = 1.0;
      jnjp1 = 1;
      jncy = 1.0;
   }

   if (zbdry == WALL) {
      k0km1 = 0;
      k0cz = 2.0;
      knkp1 = nzm1;
      kncz = 2.0;
   } else if (zbdry == PERIODIC) {
      k0km1 = nzm1 - 1;
      k0cz = 1.0;
      knkp1 = 1;
      kncz = 1.0;
   }

   for (i=0; i<nx; i++) {
      ip1 = i+1;
      im1 = i-1;
      cx = 1.0;
      if (i == 0) {
         ip1 = 1;
         im1 = i0im1;
         cx = i0cx;
      } else if (i == nxm1) {
         ip1 = inip1;
         im1 = nxm1 - 1;
         cx = incx;
      }

   for (j=0; j<ny; j++) {
      jp1 = j+1;
      jm1 = j-1;
      cy = 1.0;
      if (j == 0) {
         jp1 = 1;
         jm1 = j0jm1;
         cy = j0cy;
      } else if (j == nym1) {
         jp1 = jnjp1;
         jm1 = nym1 - 1;
         cy = jncy;
      }

   for (k=0; k<nz; k++) {
      kp1 = k+1;
      km1 = k-1;
      cz = 1.0;
      if (k == 0) {
         kp1 = 1;
         km1 = k0km1;
         cz = k0cz;
      } else if (k == nzm1) {
         kp1 = knkp1;
         km1 = nzm1 - 1;
         cz = kncz;
      }

      x_out[i][j][k] = (z_in[i][jp1][k]-z_in[i][jm1][k])*hyi*cy
                     -(y_in[i][j][kp1]-y_in[i][j][km1])*hzi*cz;

      y_out[i][j][k] = (x_in[i][j][kp1]-x_in[i][j][km1])*hzi*cz
                     -(z_in[ip1][j][k]-z_in[im1][j][k])*hxi*cx;

      z_out[i][j][k] = (y_in[ip1][j][k]-y_in[im1][j][k])*hxi*cx
                     -(x_in[i][jp1][k]-x_in[i][jm1][k])*hyi*cy;

   }
   }
   }

   /* and that's it, return 0 is all went well */
   return(0);
}


/*
 * compute curl of the input, put it in the output arrays
 */
float find_energy_3d(int nx,int ny,int nz,int xbdry,int ybdry,int zbdry,
   float ***x_in,float ***y_in,float ***z_in) {

   int nxm1 = nx-1;
   int nym1 = ny-1;
   int nzm1 = nz-1;
   int i,j,k;
   float hxi,cx,cy,cz;
   float sum = 0.;

   // one over cell size, since xsize=1, so this is always nxm1
   hxi = 1./pow((float)nxm1,3);

   for (i=0; i<nx; i++) {
      cx = 1.0;
      if (i == 0 || i == nxm1) {
         cx = 0.5;
      }

   for (j=0; j<ny; j++) {
      cy = 1.0;
      if (j == 0 || j == nym1) {
         cy = 0.5;
      }

   for (k=0; k<nz; k++) {
      cz = 1.0;
      if (k == 0 || k == nzm1) {
         cz = 0.5;
      }

      sum += (x_in[i][j][k]*x_in[i][j][k] +
              y_in[i][j][k]*y_in[i][j][k] +
              z_in[i][j][k]*z_in[i][j][k])*cx*cy*cz*hxi;

   }
   }
   }

   /* and that's it, return 0 is all went well */
   return(0.5*sum);
}


/*
 * basic 2d interpolation scheme (cloud-in-cell, or area-weighting)
 */
float interpolate_using_area_2d(int nx,int ny,int xbdry,int ybdry,float **a,float px,float py) {

   // really, this is the M2 scheme, or the cloud-in-cell scheme of
   //    bilinear interpolation

   int i,j,ii,ji,ir,jr;
   int si[2];           // start index
   FLOAT dx;
   FLOAT wt[2][2];
   FLOAT xfactor,yfactor;
   FLOAT mf;
   FLOAT out = 0.0;

   // compute the lowest corner index of the cells that the scheme needs
   // these may be negative
   si[0] = (int)(10+(nx-1)*(px)) - 10;
   // si[1] = (int)(10+(ny-1)*(py)) - 10;
   si[1] = (int)(10+(nx-1)*(py)) - 10;
   // fprintf(stdout,"location is %g %g, start index is %d %d\n",px,py,si[0],si[1]);

   // make the list of M4 weights
   dx = fabs( (px)*(nx-1) - si[0] );
   wt[0][0] = 1.0-dx;
   dx = dx-1.0;
   wt[0][1] = 1.0+dx;

   // dx = fabs( (py)*(ny-1) - si[1] );
   dx = fabs( (py)*(nx-1) - si[1] );
   wt[1][0] = 1.0-dx;
   dx = dx-1.0;
   wt[1][1] = 1.0+dx;

   // fprintf(stdout,"  weights are %g %g and %g %g\n",wt[0][0],wt[0][1],wt[1][0],wt[1][1]);

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

      /* find the wt-factor for this cell's contribution */
      mf = wt[0][ii]*wt[1][ji];

      /* apply them to the grid node in question */
      out += a[ir][jr]*mf*xfactor*yfactor;

   }}

   // fprintf(stdout,"      out is %g\n",out);

   /* all's well that ends well */
   return(out);
}


/*
 * Interpolate velocity FROM the grid using the CIC method, but in 2
 * dimensions only, and with a vector-valued field
 *
 * loc is the point at which the vector is to be evaluated (still 3d)
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
 * loc is the point at which the vector is to be evaluated (still 3d)
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
 * loc is the point at which the vector is to be evaluated (still 3d)
 * zeta makes the 2d input field of vectors
 * out is the vector value interpolated from the field
 */
float interpolate_array_using_M4p_2d(int nx,int ny,int xbdry,int ybdry,
   float ***zeta,float px,float py,int numout,float out[MAX_SCALARS]) {

   int i,j,k,ii,ji,ir,jr;
   int si[2];           // start index
   float dx;
   float m4[2][4];
   float xfactor,yfactor;
   float mf;

   for (k=0; k<numout; k++) out[k] = 0.0;

   /* compute the lowest corner index of the cells that the scheme needs */
   // these may be negative
   si[0] = (int)(10+(nx-1)*(px)) - 11;
   // si[1] = (int)(10+(ny-1)*(py)) - 11;
   si[1] = (int)(10+(nx-1)*(py)) - 11;
   // fprintf(stdout,"location is %g %g, start index is %d %d\n",px,py,si[0],si[1]);

   // make the list of M4 weights
   dx = fabs( (px)*(nx-1) - si[0] );
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
   dx = fabs( (py)*(nx-1) - si[1] );
   m4[1][0] = 0.5*pow(2-dx,2)*(1-dx);
   dx = dx-1.0;
   //m4[1][1] = 1.0 - 2.5*dx*dx + 1.5*dx*dx*dx;
   m4[1][1] = 1.0 - dx*dx*(2.5 - 1.5*dx);
   dx = fabs(dx-1.0);
   //m4[1][2] = 1.0 - 2.5*dx*dx + 1.5*dx*dx*dx;
   m4[1][2] = 1.0 - dx*dx*(2.5 - 1.5*dx);
   dx = dx+1.0;
   m4[1][3] = 0.5*pow(2-dx,2)*(1-dx);

   // fprintf(stdout,"  weights are %g %g %g %g and %g %g %g %g\n",m4[0][0],m4[0][1],m4[0][2],m4[0][3],m4[1][0],m4[1][1],m4[1][2],m4[1][3]);

   // ii,ji are the counters, i,j are the cell indexes (possibly negative)
   //   and ir,jr are the actual cells we take values from
   for (ii=0; ii<4; ii++) {
   for (ji=0; ji<4; ji++) {

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
      mf = m4[0][ii]*m4[1][ji]*xfactor*yfactor;

      /* apply them to the grid node in question */
      for (k=0; k<numout; k++) out[k] += zeta[k][ir][jr]*mf;

   }}

   /* all's well that ends well */
   return(0);
}


/*
 * Interpolate velocity FROM the grid using the CIC method
 * with a vector-valued field
 *
 * px,py,pz is the point at which the vector is to be evaluated
 * zeta is the 3d input field of vectors
 * out is the vector value interpolated from the field
 */
float interpolate_array_using_CIC_3d(int nx,int ny,int nz,
         int xbdry,int ybdry,int zbdry,float ****zeta,
         float px,float py,float pz,
         int numout,float out[MAX_SCALARS]) {

   int i,l,ci[3],oci[3];
   double x,y,z,omx,omy,omz;
   double epsilon = 1.e-5;

   // first, make sure that the point is in bounds!
   if (xbdry == PERIODIC) {
      if (px < 0.) px += 1.;
      else if (px > 1.) px -= 1.;
   } else {
      if (px < 0.) px = epsilon;
      else if (px > 1.) px = 1.-epsilon;
   }
   if (ybdry == PERIODIC) {
      if (py < 0.) py += 1.;
      else if (py > 1.) py -= 1.;
   } else {
      if (py < 0.) py = epsilon;
      else if (py > 1.) py = 1.-epsilon;
   }
   if (zbdry == PERIODIC) {
      if (pz < 0.) pz += 1.;
      else if (pz > 1.) pz -= 1.;
   } else {
      if (pz < 0.) pz = epsilon;
      else if (pz > 1.) pz = 1.-epsilon;
   }

   /* compute the lower corner index of the cell that the point xt,yt,zt is in */
   ci[0] = (int)(10 + px*(nx-1)) - 10;
   ci[1] = (int)(10 + py*(ny-1)) - 10;
   ci[2] = (int)(10 + pz*(nz-1)) - 10;

   // set "other cell index" the upper corner index
   for (i=0; i<3; i++) oci[i] = ci[i]+1;
   if (xbdry == PERIODIC && oci[0] == nx) oci[0] = 0;
   if (ybdry == PERIODIC && oci[1] == ny) oci[1] = 0;
   if (zbdry == PERIODIC && oci[2] == nz) oci[2] = 0;

   /* compute the relative location within that cell of the same point */
   x = ( px - ci[0]/(double)(nx-1) ) * (nx-1);
   y = ( py - ci[1]/(double)(ny-1) ) * (ny-1);
   z = ( pz - ci[2]/(double)(nz-1) ) * (nz-1);
   omx = 1.0-x;
   omy = 1.0-y;
   omz = 1.0-z;

   // apply them to the grid node in question
   for (l=0; l<numout; l++) {
      out[l] = omx*omy*omz*zeta[l][ci[0]][ci[1]][ci[2]] +
               omx*omy*z*  zeta[l][ci[0]][ci[1]][oci[2]] +
               omx*y*omz*  zeta[l][ci[0]][oci[1]][ci[2]] +
               omx*y*  z*  zeta[l][ci[0]][oci[1]][oci[2]] +
               x*omy*omz*  zeta[l][oci[0]][ci[1]][ci[2]] +
               x*omy*  z*  zeta[l][oci[0]][ci[1]][oci[2]] +
               x*  y*omz*  zeta[l][oci[0]][oci[1]][ci[2]] +
               x*  y*  z*  zeta[l][oci[0]][oci[1]][oci[2]];
   }

   // done
   return(0);
}


/*
 * Interpolate velocity FROM the grid using the M4' method, in 3
 * dimensions, and with a vector-valued field
 *
 * loc is the point at which the vector is to be evaluated (still 3d)
 * zeta makes the 2d input field of vectors
 * out is the vector value interpolated from the field
 */
float interpolate_array_using_M4p_3d(int nx,int ny,int nz,
         int xbdry,int ybdry,int zbdry,float ****zeta,
         float px,float py,float pz,int numout,float out[MAX_SCALARS]) {

   int i,j,k,l,ii,ji,ki,ir,jr,kr;
   int si[3];           // start index
   float dx;
   float m4[3][4];
   float xfactor,yfactor,zfactor;
   float mf;

   for (k=0; k<numout; k++) out[k] = 0.0;

   /* compute the lowest corner index of the cells that the scheme needs */
   // these may be negative
   si[0] = (int)(10+(nx-1)*(px)) - 11;
   si[1] = (int)(10+(ny-1)*(py)) - 11;
   si[2] = (int)(10+(nz-1)*(pz)) - 11;
   // fprintf(stdout,"location is %g %g %g, start index is %d %d %d\n",px,py,pz,si[0],si[1],si[2]);
   // fflush(stdout);

   // make the list of M4 weights
   dx = fabs( (px)*(nx-1) - si[0] );
   m4[0][0] = 0.5*pow(2-dx,2)*(1-dx);
   dx = dx-1.0;
   m4[0][1] = 1.0 - 2.5*dx*dx + 1.5*dx*dx*dx;
   dx = fabs(dx-1.0);
   m4[0][2] = 1.0 - 2.5*dx*dx + 1.5*dx*dx*dx;
   dx = dx+1.0;
   m4[0][3] = 0.5*pow(2-dx,2)*(1-dx);

   // dx = fabs( (py)*(ny-1) - si[1] );
   dx = fabs( (py)*(nx-1) - si[1] );
   m4[1][0] = 0.5*pow(2-dx,2)*(1-dx);
   dx = dx-1.0;
   m4[1][1] = 1.0 - 2.5*dx*dx + 1.5*dx*dx*dx;
   dx = fabs(dx-1.0);
   m4[1][2] = 1.0 - 2.5*dx*dx + 1.5*dx*dx*dx;
   dx = dx+1.0;
   m4[1][3] = 0.5*pow(2-dx,2)*(1-dx);

   dx = fabs( (pz)*(nx-1) - si[2] );
   m4[2][0] = 0.5*pow(2-dx,2)*(1-dx);
   dx = dx-1.0;
   m4[2][1] = 1.0 - 2.5*dx*dx + 1.5*dx*dx*dx;
   dx = fabs(dx-1.0);
   m4[2][2] = 1.0 - 2.5*dx*dx + 1.5*dx*dx*dx;
   dx = dx+1.0;
   m4[2][3] = 0.5*pow(2-dx,2)*(1-dx);

   // fprintf(stdout,"  weights are %g %g %g %g and %g %g %g %g\n",m4[0][0],m4[0][1],m4[0][2],m4[0][3],m4[1][0],m4[1][1],m4[1][2],m4[1][3]);
   // fflush(stdout);

   // ii,ji are the counters, i,j are the cell indexes (possibly negative)
   //   and ir,jr are the actual cells we take values from
   for (ii=0; ii<4; ii++) {
   for (ji=0; ji<4; ji++) {
   for (ki=0; ki<4; ki++) {

      i = ii+si[0];
      j = ji+si[1];
      k = ki+si[2];

      /* in order to properly reflect vorticity off walls, use xyzfactors */
      xfactor = 1.0;
      yfactor = 1.0;
      zfactor = 1.0;

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

      /* near the min or max y bounds */
      if (zbdry == PERIODIC) {
         kr = (k+(nz-1))%(nz-1);
      } else {
         kr = k;
         if (k<0) {
            kr = -k;
            // yfactor *= -1.0;
         } else if (k>(nz-1)) {
            kr = 2*(nz-1) - k;
            // yfactor *= -1.0;
         }
      }

      // fprintf(stdout,"index %d %d %d, m-index %d %d %d\n",ir,jr,kr,ii,ji,ki); fflush(stdout);

      /* find the m4-factor for this cell's contribution */
      mf = m4[0][ii]*m4[1][ji]*m4[2][ki]*xfactor*yfactor*zfactor;

      /* apply them to the grid node in question */
      for (l=0; l<numout; l++) {
         // fprintf(stdout," var %d\n",l); fflush(stdout);
         out[l] += zeta[l][ir][jr][kr]*mf;
      }

   }}}

   /* all's well that ends well */
   return(0);
}


/*
 * Interpolate velocity FROM the grid using the M4' method, but in 2
 * dimensions only, and with a scalar-valued field
 *
 * loc is the point at which the vector is to be evaluated (still 3d)
 * zeta makes the 2d input field of vectors
 * out is the scalar value interpolated from the field
 */
float interpolate_using_M4p_2d(int nx,int ny,int xbdry,int ybdry,float **zeta,float px,float py) {

   int i,j,ii,ji,ir,jr;
   int si[2];           // start index
   FLOAT dx;
   FLOAT m4[2][4];
   FLOAT xfactor,yfactor;
   FLOAT mf;
   FLOAT out = 0.0;

   /* compute the lowest corner index of the cells that the scheme needs */
   // these may be negative
   si[0] = (int)(10+(nx-1)*(px)) - 11;
   // si[1] = (int)(10+(ny-1)*(py)) - 11;
   si[1] = (int)(10+(nx-1)*(py)) - 11;
   // fprintf(stdout,"location is %g %g, start index is %d %d\n",px,py,si[0],si[1]);

   // make the list of M4 weights
   dx = fabs( (px)*(nx-1) - si[0] );
   m4[0][0] = 0.5*pow(2-dx,2)*(1-dx);
   dx = dx-1.0;
   m4[0][1] = 1.0 - 2.5*dx*dx + 1.5*dx*dx*dx;
   dx = fabs(dx-1.0);
   m4[0][2] = 1.0 - 2.5*dx*dx + 1.5*dx*dx*dx;
   dx = dx+1.0;
   m4[0][3] = 0.5*pow(2-dx,2)*(1-dx);

   // dx = fabs( (py)*(ny-1) - si[1] );
   dx = fabs( (py)*(nx-1) - si[1] );
   m4[1][0] = 0.5*pow(2-dx,2)*(1-dx);
   dx = dx-1.0;
   m4[1][1] = 1.0 - 2.5*dx*dx + 1.5*dx*dx*dx;
   dx = fabs(dx-1.0);
   m4[1][2] = 1.0 - 2.5*dx*dx + 1.5*dx*dx*dx;
   dx = dx+1.0;
   m4[1][3] = 0.5*pow(2-dx,2)*(1-dx);

   // fprintf(stdout,"  weights are %g %g %g %g and %g %g %g %g\n",m4[0][0],m4[0][1],m4[0][2],m4[0][3],m4[1][0],m4[1][1],m4[1][2],m4[1][3]);

   // ii,ji are the counters, i,j are the cell indexes (possibly negative)
   //   and ir,jr are the actual cells we take values from
   for (ii=0; ii<4; ii++) {
   for (ji=0; ji<4; ji++) {

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
      mf = m4[0][ii]*m4[1][ji];

      /* apply them to the grid node in question */
      out += zeta[ir][jr]*mf*xfactor*yfactor;

   }}

   // fprintf(stdout,"      out is %g\n",out);

   /* all's well that ends well */
   return(out);
}


/*
 * Interpolate velocity FROM the grid using the CIC method, but in 2
 * dimensions only, and with a vector-valued field
 *
 * loc is the point at which the vector is to be evaluated (still 3d)
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
 * loc is the point at which the vector is to be evaluated (still 3d)
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
 * loc is the point at which the vector is to be evaluated (still 3d)
 * zeta makes the 2d input field of vectors
 * out is the vector value interpolated from the field
 */
int interpolate_vel_using_M4p_2d(int nx,int ny,int xbdry,int ybdry,float **ua,float **va,float px,float py,float *u,float *v) {

   int i,j,ii,ji,ir,jr;
   int si[2];           // start index
   FLOAT dx;
   FLOAT m4[2][4];
   FLOAT xfactor,yfactor;
   FLOAT mf;
   FLOAT out = 0.0;

   /* compute the lowest corner index of the cells that the scheme needs */
   // these may be negative
   si[0] = (int)(10+(nx-1)*(px)) - 11;
   // si[1] = (int)(10+(ny-1)*(py)) - 11;
   si[1] = (int)(10+(nx-1)*(py)) - 11;
   // fprintf(stdout,"location is %g %g, start index is %d %d\n",px,py,si[0],si[1]);

   // make the list of M4 weights
   dx = fabs( (px)*(nx-1) - si[0] );
   m4[0][0] = 0.5*pow(2-dx,2)*(1-dx);
   dx = dx-1.0;
   m4[0][1] = 1.0 - 2.5*dx*dx + 1.5*dx*dx*dx;
   dx = fabs(dx-1.0);
   m4[0][2] = 1.0 - 2.5*dx*dx + 1.5*dx*dx*dx;
   dx = dx+1.0;
   m4[0][3] = 0.5*pow(2-dx,2)*(1-dx);

   // dx = fabs( (py)*(ny-1) - si[1] );
   dx = fabs( (py)*(nx-1) - si[1] );
   m4[1][0] = 0.5*pow(2-dx,2)*(1-dx);
   dx = dx-1.0;
   m4[1][1] = 1.0 - 2.5*dx*dx + 1.5*dx*dx*dx;
   dx = fabs(dx-1.0);
   m4[1][2] = 1.0 - 2.5*dx*dx + 1.5*dx*dx*dx;
   dx = dx+1.0;
   m4[1][3] = 0.5*pow(2-dx,2)*(1-dx);

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
      mf = xfactor*yfactor*m4[0][ii]*m4[1][ji];

      /* apply them to the grid node in question */
      *(u) += ua[ir][jr]*mf;
      *(v) += va[ir][jr]*mf;

   }}

   // fprintf(stdout,"      out is %g\n",out);

   /* all's well that ends well */
   return(0);
}


/*
 * Copy an array
 */
int copy_2d_field(float **from, float **to, int nx, int ny) {
   int i,j;

   for (i=0;i<nx;i++) for (j=0;j<ny;j++) {
      to[i][j] = from[i][j];
   }

   return(0);
}


/*
 * allocate memory for a one-dimensional array of float
 */
float* allocate_1d_array_f(int nx) {
 
   float *array = (float *)malloc(nx * sizeof(float));
 
   return(array);
}


int free_1d_array_f(float* array){
   free(array);
   return(0);
}


/*
 * allocate memory for a two-dimensional array of float
 */
float** allocate_2d_array_f(int nx,int ny) {
 
   int i;
   float **array = (float **)malloc(nx * sizeof(float *));
 
   array[0] = (float *)malloc(nx * ny * sizeof(float));
   for (i=1; i<nx; i++)
      array[i] = array[0] + i * ny;
 
   return(array);
}


int free_2d_array_f(float** array){
   free(array[0]);
   free(array);
   return(0);
}


/*
 * allocate memory for a three-dimensional array of floats
 */
float*** allocate_3d_array_f(int nx, int ny, int nz) {
 
   int i,j;
   float ***array = (float ***)malloc(nx * sizeof(float **));
 
   array[0] = (float **)malloc(nx * ny * sizeof(float *));
   array[0][0] = (float *)malloc(nx * ny * nz * sizeof(float));
 
   for (i=1; i<nx; i++)
      array[i] = array[0] + i * ny;
 
   for (i=0; i<nx; i++) {
      if (i!=0)
         array[i][0] = array[0][0] + i * ny * nz;
      for (j=1; j<ny; j++)
         array[i][j] = array[i][0] + j * nz;
   }
 
   return(array);
}


int free_3d_array_f(float*** array){
   free(array[0][0]);
   free(array[0]);
   free(array);
   return(0);
}


/* ------------------------------------------------ *
 * gaussian -- generates a gaussian random variable *
 *             with mean a and standard deviation d *
 * ------------------------------------------------ */
double gaussian(double a,double d) {

   static double t = 0.0;
   double x,v1,v2,r;

   if (t == 0) {
      do {
         // v1 = 2.0 * rand() - 1.0;
         // v2 = 2.0 * rand() - 1.0;
         v1 = 2.*rand()/(RAND_MAX+1.0) - 1.;
         v2 = 2.*rand()/(RAND_MAX+1.0) - 1.;
         r = v1 * v1 + v2 * v2;
      } while (r>=1.0);
      r = sqrt((-2.0*log(r))/r);
      t = v2*r;
      return(a+v1*r*d);
   } else {
      x = t;
      t = 0.0;
      return(a+x*d);
   }
}
