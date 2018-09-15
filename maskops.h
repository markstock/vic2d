/*
 * VIC-MOC - maskops.h - operations on dynamics masks
 *
 * Copyright 2014 Mark J. Stock mstock@umich.edu
 */

#pragma once

typedef struct block_record {
   float tstart;
   float tlen;
   int startx,endx,starty,endy;
   float r,g,b;
   int mask;
} BLOCK;

#define MAXBLOCKS 10000
BLOCK block[MAXBLOCKS];
int nblocks = 0;

void populate_block_array (int, int);
void update_mask_with_blocks_1 (float**, float***, int, int, float);
void update_mask_with_blocks_2 (float**, float***, int, int, float, float);
void overlay_mask (int, int, float**);

