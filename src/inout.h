/*
 * VIC-MOC - inout.h - routines for input and output
 *
 * Copyright 2004-7 Mark J. Stock mstock@umich.edu
 */

#pragma once

int write_output(char*,int,int,float**,float,float,int);
int write_png (char*, int, int, int, int, float**, float, float, float**, float, float, float**, float, float);
int read_png_res (char *infile, int *hgt, int *wdt);
int read_png (char*, int, int, int, int, float, int, float**, float, float, float**, float, float, float**, float, float);
int write_output_3d(char*,int,int,int,float***,float,float,int,int);
int write_3d_vtk(char*,int,int,int,float***,float***,float***);
int write_output_particles_rad(char*,int,float**,float**,float*);


