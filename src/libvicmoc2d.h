/*
 * VIC-MOC - libvicmoc2d.h - structures and defines for 2d lib
 *
 * Copyright 2004-7,20,23 Mark J. Stock <markjstock@gmail.com>
 */

int find_vels_2d (int,int,int,int,int,int,int,float*,float*,float**,float**,float**,const int,float**,const float);
int add_sharp_circular_blob(int,int,int,int,float**,float,float,float,float);
int add_smooth_circular_blob(int,int,int,int,float**,float,float,float,float);
int find_shear_magnitude (int, int, int, int, float**, float, float**, float, float **);
float step_forward_2d(int,int,int,int,int,int,int,int,float*,float*,int,int,float***,float***,float***,const int,float**,const float,float*,int,float*,float,int,float,float,float***);

