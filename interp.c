/*
 * VIC-MOC - advection and interpolation routine for vic2d/vic3d
 *
 * Copyright 2011 Mark J. Stock mstock@umich.edu
 */

#include "vicmoc.h"
#ifdef _OPENMP
#include <omp.h>
#endif

// DIM is dimension, 2 or 3
// SUPER 1 means no subsampling, 2 means 2x2 subsampling
// METH is CIC, TSC, M4P or M6P
// ORDER is 1, 2, or 4


#if DIM==2
int moc_advect_DIMd_SUPER_METH_ORDERth (int,int,int,int,float**,float**,float***,float***,float);
#else
int moc_advect_DIMd_SUPER_METH_ORDERth (int,int,int,int,int,int,float***,float***,float***,float***,float***,float);
#endif


/*
 * moc_advect_2d is the implicit scheme for updating the vorticity values
 */
int moc_advect_DIMd_SUPER_METH_ORDERth
                  (int nx,int ny,
#                  if DIM==3
                      int nz,
#                  endif
                   int xbdry,int ybdry,
#                  if DIM==3
                      int zbdry,
#                  endif
#                  if DIM==2
                     float **u,float **v,
#                  else
                     float ***u,float ***v,float ***w,
#                  endif
                   float ***in,float ***out,
                   float dt) {

   int i,j,k;
   int sc_cnt;
   float xf,yf,accx,accy,px,py,newx,newy;
#if DIM==3
   float zf,accz,pz,newz;
#endif
   float u0,u1,u2,u3,v0,v1,v2,v3;
   float dx,oodt,temp;
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
      //fprintf(stderr,"%d  %ld  %ld\n",i,in[i],out[i]);
      //fflush(stderr);
      if (in[i] != NULL) {
         tempin[sc_cnt] = in[i];
         tempout[sc_cnt] = out[i];
         //fprintf(stderr,"   %ld  %ld\n",tempin[sc_cnt],tempout[sc_cnt]);
         //fflush(stderr);
         sc_cnt++;
      }
   }
   //for (i=0; i<MAX_SCALARS; i++) {
   //   fprintf(stderr,"%ld  %ld\n",in[i],tempin[i]);
   //}
   //fflush(stderr);

   // fprintf(stderr,"in moc_advect_2d\n"); fflush(stderr);
#if DIM==2
   if (nx > ny) {
      xf = 1.0;
      yf = (float)(ny-1)/(float)(nx-1);
      dx = 1.0 / (float)(nx-1);
   } else {
      xf = (float)(nx-1)/(float)(ny-1);
      yf = 1.0;
      dx = 1.0 / (float)(ny-1);
   }
#else
   if (nx > ny) {
      if (nx > nz) {
         xf = 1.0;
         yf = (float)(ny-1)/(float)(nx-1);
         zf = (float)(nz-1)/(float)(nx-1);
         dx = 1.0 / (float)(nx-1);
      } else {
         xf = (float)(nx-1)/(float)(nz-1);
         yf = (float)(ny-1)/(float)(nz-1);
         zf = 1.0;
         dx = 1.0 / (float)(nz-1);
      }
   } else {
      if (ny > nz) {
         xf = (float)(nx-1)/(float)(ny-1);
         yf = 1.0;
         zf = (float)(nz-1)/(float)(ny-1);
         dx = 1.0 / (float)(ny-1);
      } else {
         xf = (float)(nx-1)/(float)(nz-1);
         yf = (float)(ny-1)/(float)(nz-1);
         zf = 1.0;
         dx = 1.0 / (float)(nz-1);
      }
   }
#endif

   // for each point, march backwards in space and time to find the value
#  if DIM==2
      #pragma omp parallel for private(i,px,j,py,newx,newy,k,outvals)
#  else
      #pragma omp parallel for private(i,px,j,py,k,pz,newx,newy,newz,outvals)
#  endif

   for (i=0; i<nx; i++) {
   px = (float)i * dx;
#  if SUPER>1 // subsample
   for (ii=0; ii<SUPER; ii++) {
   pxx = px + ((ii+0.5)/SUPER - 0.5) * dx;
#  else
   pxx = px;
#  endif

    for (j=0; j<ny; j++) {
    py = (float)j * dx;
#   if SUPER>1 // subsample
    for (jj=0; jj<SUPER; jj++) {
    pyy = py + ((jj+0.5)/SUPER - 0.5) * dx;
#   else
    pyy = py;
#   endif

#    if DIM==3
     for (k=0; k<nz; k++) {
     pz = (float)k * dx;
#    if SUPER>1 // subsample
     for (kk=0; kk<SUPER; kk++) {
     pzz = pz + ((kk+0.5)/SUPER - 0.5) * dx;
#    else
     pzz = pz;
#    endif
#    endif


#     if SUPER>1 // subsample
#      if DIM==2
        // find vel at px,py
        interpolate_vel_using_METH_2d(nx,ny,xbdry,ybdry,u,v,pxx,pyy,&u0,&v0);
#      else
        // find vel at px,py,pz
        interpolate_vel_using_METH_3d(nx,ny,nz,xbdry,ybdry,zbdry,u,v,w,pxx,pyy,pzz,&u0,&v0,&w0);
#      endif
#     else
#      if DIM==2
        // find vel at px,py
        u0 = u[i][j];
        v0 = v[i][j];
#      else
        // find vel at px,py,pz
        u0 = u[i][j][k];
        v0 = v[i][j][k];
        w0 = w[i][j][k];
#      endif
#     endif

#     if DIM==2
         newx = px-dt*u[i][j];
         newy = py-dt*v[i][j];
#     else
         newx = px-dt*u[i][j][k];
         newy = py-dt*v[i][j][k];
         newz = pz-dt*w[i][j][k];
#     endif

         if (xbdry == WALL || xbdry == OPEN) {
            if (newx > xf-EPSILON) newx = xf-EPSILON;
            if (newx < 0.0+EPSILON) newx = EPSILON;
         }
         if (ybdry == WALL || ybdry == OPEN) {
            if (newy > yf-EPSILON) newy = yf-EPSILON;
            if (newy < 0.0+EPSILON) newy = EPSILON;
         }
#     if DIM==3
         if (zbdry == WALL || zbdry == OPEN) {
            if (newz > zf-EPSILON) newz = zf-EPSILON;
            if (newz < 0.0+EPSILON) newz = EPSILON;
         }
#     endif

#     if METH==CIC
#     elif METH==TSC
#     elif METH==M4P
#     else
#     endif



         // last thing---interpolate all other variables
#        if DIM==2
            for (l=0; l<sc_cnt; l++) tempout[l][i][j] = outvals[l];
#        else
            for (l=0; l<sc_cnt; l++) tempout[l][i][j][k] = outvals[l];
#        endif

#        if DIM==3
         }
#        endif



#   if DIM==3
     }
#    if SUPER>1 // subsample
     }
#    endif
#   endif

#   if SUPER>1 // subsample
    }
#   endif
    }

#  if SUPER>1 // subsample
   }
#  endif
   }


// second order

            // find vel at newx,newy
            if (interp == cic)
               interpolate_vel_using_CIC_2d(nx,ny,xbdry,ybdry,u,v,newx,newy,&u1,&v1);
            else if (interp == tsc)
               interpolate_vel_using_TSC_2d(nx,ny,xbdry,ybdry,u,v,newx,newy,&u1,&v1);
            else
               interpolate_vel_using_M4p_2d(nx,ny,xbdry,ybdry,u,v,newx,newy,&u1,&v1);

#ifdef ORDER>1
            // find position back 1 step using average of two velocities
            newx = px-dt*0.5*(u0+u1);
            if (xbdry == WALL || xbdry == OPEN) {
               if (newx > xf-EPSILON) newx = xf-EPSILON;
               if (newx < 0.0+EPSILON) newx = EPSILON;
            }
            newy = py-dt*0.5*(v0+v1);
            if (ybdry == WALL || ybdry == OPEN) {
               if (newy > yf-EPSILON) newy = yf-EPSILON;
               if (newy < 0.0+EPSILON) newy = EPSILON;
            }

            // are we able to calculate the acceleration here?

            // finally, return scalar-valued vorticity from there
            if (interp == cic)
               interpolate_array_using_CIC_2d(nx,ny,xbdry,ybdry,tempin,newx,newy,sc_cnt,outvals);
            else if (interp == tsc)
               interpolate_array_using_TSC_2d(nx,ny,xbdry,ybdry,tempin,newx,newy,sc_cnt,outvals);
            else
               interpolate_array_using_M4p_2d(nx,ny,xbdry,ybdry,tempin,newx,newy,sc_cnt,outvals);
            for (k=0; k<sc_cnt; k++) tempout[k][i][j] = outvals[k];
         }
      }
   } else if (order == 4) {
      // use a four-step method, RK4
      #pragma omp parallel for private(i,px,j,py,u0,u1,u2,u3,v0,v1,v2,v3,newx,newy,accx,accy,k,outvals)
      for (i=0;i<nx;i++) {
         px = (float)i * dx;
         for (j=0;j<ny;j++) {

            py = (float)j * dx;
            // find vel at px,py (k_1)
            u0 = u[i][j];
            v0 = v[i][j];

            // find position 1 explicit Euler step backwards
            newx = px-dt*0.5*u0;
            if (xbdry == WALL || xbdry == OPEN) {
               if (newx > xf-EPSILON) newx = xf-EPSILON;
               if (newx < 0.0+EPSILON) newx = EPSILON;
            }
            newy = py-dt*0.5*v0;
            if (ybdry == WALL || ybdry == OPEN) {
               if (newy > yf-EPSILON) newy = yf-EPSILON;
               if (newy < 0.0+EPSILON) newy = EPSILON;
            }
            // find vel at newx,newy (k_2)
            if (interp == cic)
               interpolate_vel_using_CIC_2d(nx,ny,xbdry,ybdry,u,v,newx,newy,&u1,&v1);
            else if (interp == tsc)
               interpolate_vel_using_TSC_2d(nx,ny,xbdry,ybdry,u,v,newx,newy,&u1,&v1);
            else
               interpolate_vel_using_M4p_2d(nx,ny,xbdry,ybdry,u,v,newx,newy,&u1,&v1);

            // find position 2 explicit Euler step backwards
            newx = px-dt*0.5*u1;
            if (xbdry == WALL || xbdry == OPEN) {
               if (newx > xf-EPSILON) newx = xf-EPSILON;
               if (newx < 0.0+EPSILON) newx = EPSILON;
            }
            newy = py-dt*0.5*v1;
            if (ybdry == WALL || ybdry == OPEN) {
               if (newy > yf-EPSILON) newy = yf-EPSILON;
               if (newy < 0.0+EPSILON) newy = EPSILON;
            }
            // find vel at newx,newy (k_3)
            if (interp == cic)
               interpolate_vel_using_CIC_2d(nx,ny,xbdry,ybdry,u,v,newx,newy,&u2,&v2);
            else if (interp == tsc)
               interpolate_vel_using_TSC_2d(nx,ny,xbdry,ybdry,u,v,newx,newy,&u2,&v2);
            else
               interpolate_vel_using_M4p_2d(nx,ny,xbdry,ybdry,u,v,newx,newy,&u2,&v2);

            // find position 3 explicit Euler step backwards
            newx = px-dt*u2;
            if (xbdry == WALL || xbdry == OPEN) {
               if (newx > xf-EPSILON) newx = xf-EPSILON;
               if (newx < 0.0+EPSILON) newx = EPSILON;
            }
            newy = py-dt*v2;
            if (ybdry == WALL || ybdry == OPEN) {
               if (newy > yf-EPSILON) newy = yf-EPSILON;
               if (newy < 0.0+EPSILON) newy = EPSILON;
            }
            // find vel at newx,newy (k_4)
            if (interp == cic)
               interpolate_vel_using_CIC_2d(nx,ny,xbdry,ybdry,u,v,newx,newy,&u3,&v3);
            else if (interp == tsc)
               interpolate_vel_using_TSC_2d(nx,ny,xbdry,ybdry,u,v,newx,newy,&u3,&v3);
            else
               interpolate_vel_using_M4p_2d(nx,ny,xbdry,ybdry,u,v,newx,newy,&u3,&v3);

            // find position back 1 step using average of four velocities
            accx = 0.16666667*(u0+u3+2.*(u1+u2));
            newx = px-dt*accx;
            if (xbdry == WALL || xbdry == OPEN) {
               if (newx > xf-EPSILON) newx = xf-EPSILON;
               if (newx < 0.0+EPSILON) newx = EPSILON;
            }
            accy = 0.16666667*(v0+v3+2.*(v1+v2));
            newy = py-dt*accy;
            if (ybdry == WALL || ybdry == OPEN) {
               if (newy > yf-EPSILON) newy = yf-EPSILON;
               if (newy < 0.0+EPSILON) newy = EPSILON;
            }

            // gather an estimate of the acceleration (u1-u0)/dt?
            accx = (u0-accx)*oodt;
            accy = (v0-accy)*oodt;
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
               interpolate_array_using_M4p_2d(nx,ny,xbdry,ybdry,tempin,newx,newy,sc_cnt,outvals);
            for (k=0; k<sc_cnt; k++) tempout[k][i][j] = outvals[k];
         }
      }
   } else {
      fprintf(stderr,"Method-of-characteristics advection order %d unsupported.\n",order);
      fprintf(stderr," Try 1 or 2.\n");
      fprintf(stderr," Quitting.\n");
      exit(0);
   }

   return(0);
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
   float ***zeta,float px,float py,int numout,float out[MAX_SCALARS]) {

   int i,j,k,ii,ji,ir,jr,nm1;
   int si[2];           // start index
   float dx;
   float m4[2][4];
   float xfactor,yfactor;
   float mf;

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
   // fprintf(stdout,"location is %g %g, start index is %d %d\n",px,py,si[0],si[1]);

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
 * Interpolate velocity FROM the grid using the M4' method, but in 2
 * dimensions only, and with a scalar-valued field
 *
 * loc is the point at which the vector is to be evaluated
 * zeta makes the 2d input field of vectors
 * out is the scalar value interpolated from the field
 */
float interpolate_using_M4p_2d(int nx,int ny,int xbdry,int ybdry,float **zeta,float px,float py) {

   int i,j,ii,ji,ir,jr,nm1;
   int si[2];           // start index
   FLOAT dx;
   FLOAT m4[2][4];
   FLOAT xfactor,yfactor;
   FLOAT mf;
   FLOAT out = 0.0;

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
   dx = fabs( (px)*(nm1) - si[0] );
   m4[0][0] = 0.5*pow(2-dx,2)*(1-dx);
   dx = dx-1.0;
   m4[0][1] = 1.0 - 2.5*dx*dx + 1.5*dx*dx*dx;
   dx = fabs(dx-1.0);
   m4[0][2] = 1.0 - 2.5*dx*dx + 1.5*dx*dx*dx;
   dx = dx+1.0;
   m4[0][3] = 0.5*pow(2-dx,2)*(1-dx);

   dx = fabs( (py)*(nm1) - si[1] );
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
int interpolate_vel_using_M4p_2d(int nx,int ny,int xbdry,int ybdry,float **ua,float **va,float px,float py,float *u,float *v) {

   int i,j,ii,ji,ir,jr,nm1;
   int si[2];           // start index
   FLOAT dx;
   FLOAT m4[2][4];
   FLOAT xfactor,yfactor;
   FLOAT mf;
   FLOAT out = 0.0;

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
      //mf = xfactor*yfactor*m4[0][ii]*m4[1][ji];
      mf = m4[0][ii]*m4[1][ji];

      /* apply them to the grid node in question */
      *(u) += ua[ir][jr]*mf;
      *(v) += va[ir][jr]*mf;

   }}

   // fprintf(stdout,"      out is %g\n",out);

   /* all's well that ends well */
   return(0);
}


