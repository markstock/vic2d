/*
 * VIC-MOC - maskops.h - operations on dynamic masks
 *
 * Copyright 2014,20 Mark J. Stock <mstock@umich.edu>
 */

#pragma once


void overlay_mask (int, int, float**);
void set_mask_from_temporal (int, int, int, float**, char*, int, float, int, float, float, float);

void mod_mask_with_vel (int, int, float, float**, float**, float**);
void mod_mask_with_vort (int, int, float, float**, float**);

void populate_block_array (int, int);
void update_mask_with_blocks_1 (float**, float***, int, int, float);
void update_mask_with_blocks_2 (float**, float***, int, int, float, float);

