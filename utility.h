/*
 * VIC-MOC - utility.h - the library for the 2D routines
 *
 * Copyright 2004-10 Mark J. Stock mstock@umich.edu
 */

#pragma once

// allocation and utility routines
int copy_2d_field (float**,float**,int,int);
float* allocate_1d_array_f (long int);
int free_1d_array_f (float*);
float** allocate_2d_array_f (int,int);
int free_2d_array_f (float**);
float*** allocate_3d_array_f(int,int,int);
int free_3d_array_f(float***);
float** flatten_to_2d (float***, int, int, int, int);
int** allocate_2d_array_i (int,int);
int free_2d_array_i (int**);
double gaussian(double, double);
