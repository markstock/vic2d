/*
 * VIC-MOC - utility.c - the library for the 2D routines
 *
 * Copyright 2004-10 Mark J. Stock <mstock@umich.edu>
 */

#include <stdlib.h>
#include <stdio.h>
//#ifdef __LINUX__
#include <malloc.h>
//#endif
#include <math.h>

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
float* allocate_1d_array_f(long int nx) {
 
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
 * take a 3D array and integrate it to make a 2D array along axis
 */

float** flatten_to_2d (float ***a, int axis, int nx, int ny, int nz) {

   int i,j,k;
   float **retArry;

   if (axis == 0) {
      // integrate over x axis
      retArry = allocate_2d_array_f (ny,nz);
         for (j=0;j<ny;j++) {
            for (k=0;k<nz;k++) {
               retArry[j][k] = 0.;
            }
         }
      for (i=0;i<nx;i++) {
         for (j=0;j<ny;j++) {
            for (k=0;k<nz;k++) {
               retArry[j][k] += a[i][j][k];
            }
         }
      }
         for (j=0;j<ny;j++) {
            for (k=0;k<nz;k++) {
               retArry[j][k] /= (float)nx;
            }
         }

   } else if (axis == 1) {
      // integrate over y axis
      retArry = allocate_2d_array_f (nz,nx);
      for (i=0;i<nx;i++) {
            for (k=0;k<nz;k++) {
               retArry[k][i] = 0.;
            }
      }
      for (i=0;i<nx;i++) {
         for (j=0;j<ny;j++) {
            for (k=0;k<nz;k++) {
               retArry[k][i] += a[i][j][k];
            }
         }
      }
      for (i=0;i<nx;i++) {
            for (k=0;k<nz;k++) {
               retArry[k][i] /= (float)ny;
            }
      }
   } else if (axis == 2) {
      // integrate over z axis
      retArry = allocate_2d_array_f (nx,ny);
      for (i=0;i<nx;i++) {
         for (j=0;j<ny;j++) {
               retArry[i][j] = 0.;
         }
      }
      for (i=0;i<nx;i++) {
         for (j=0;j<ny;j++) {
            for (k=0;k<nz;k++) {
               retArry[i][j] += a[i][j][k];
            }
         }
      }
      for (i=0;i<nx;i++) {
         for (j=0;j<ny;j++) {
               retArry[i][j] /= (float)nz;
         }
      }
   } else {
      fprintf(stderr,"ERROR (flatten_to_2d): invalid axis (%d)",axis);
      exit(1);
   }

   return(retArry);
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


/*
 * allocate memory for a two-dimensional array of ints
 */
int** allocate_2d_array_i(int nx,int ny) {
 
   int i;
   int **array = (int **)malloc(nx * sizeof(int *));
 
   array[0] = (int *)malloc(nx * ny * sizeof(int));
   for (i=1; i<nx; i++)
      array[i] = array[0] + i * ny;
 
   return(array);
}


int free_2d_array_i(int** array){
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

