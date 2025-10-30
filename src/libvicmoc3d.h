/*
 * VIC-MOC - libvicmoc3d.h - structures and defines for 3d lib
 *
 * Copyright 2004-7,20,23 Mark J. Stock <markjstock@gmail.com>
 */

float find_vels_3d(int,int,int,int,int,int,int,float***,float***,float***,float***,float***,float***);
int add_smooth_spherical_blob(int,int,int,int,int,int,float***,float,float,float,float,float);
int add_singular_blob_3d(int,int,int,int,int,int,float***,float,float,float,float,float);
int add_cube_3d(int,int,int,int,int,int,float***,float,float,float,float,float);
int add_vortex_ring_3d(int,int,int,int,int,int,float***,float***,float***,float,float,float,float,float,float,float,float,float);
int add_smooth_circular_blob_3d(int,int,int,int,int,int,int,float***,float,float,float,float);
int explicit_particle_move_3d(int,int,int,int,int,int,float****,float,float,int,float**,float**);
float step_forward_3d(int,int,int,int,int,int,int,int,float****,float****,float****,int,int*,float*,float,float);
int make_solenoidal_3d(int,int,int,int,int,int,float***,float***,float***);
float find_energy_3d(int,int,int,int,int,int,float***,float***,float***);
float find_vmax(float***, float***, float***,int,int,int);
