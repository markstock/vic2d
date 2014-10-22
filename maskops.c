#include "vicmoc.h"

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

//
// create routine to pre-generate and save a long list of "add" and "subtract" blocks
//    and make sure to force each new block to uncover/cover something new
// incorporate some routine to grab colors from a still frame (have that already?)
// generate or find some 
//

void populate_block_array (int nx, int ny) {

   // first, prepare for colors

   // the color image from which to grab colors
   char colorsrcfilename[255];
   sprintf(colorsrcfilename,"color_source.png");
   const int cnx = 300;
   const int cny = 169;

   // allocate space and read image
   float **c[3];
   c[0] = allocate_2d_array_f(cnx,cny);
   c[1] = allocate_2d_array_f(cnx,cny);
   c[2] = allocate_2d_array_f(cnx,cny);
   read_png (colorsrcfilename,cnx,cny,TRUE,FALSE,1.0,FALSE,
             c[0],0.0,1.0,c[1],0.0,1.0,c[2],0.0,1.0);

   // later on, grab colors with
   //(void) get_random_color (c,cnx,cny,thisc);

   // allocate space for mask queries
   int **masked;
   masked = allocate_2d_array_i(nx,ny);
   for (int ix=0; ix<nx; ix++) for (int iy=0; iy<ny; iy++) masked[ix][iy] = TRUE;

   // then, generate the blocks
   int icnt = 0;

   // time center and half-width of sequence of blocks that will open the domain
   const float topencenter = 55.;
   const float topenwide = 60.;
   const float topenmag = 9.5;
   // time center and half-width of sequence of blocks that will mask the domain
   const float tmaskcenter = 170.;
   const float tmaskwide = 45.;
   const float tmaskmag = 13.5;
   const int stop_time = 310;
   // dimensions of blocks
   const float sizebase = (float)nx*30./2049.;
   const float sizerand = (float)nx*100./2049.;

   // blocks that open up the space (mostly)
   for (int i=(int)(MAX(0.,topencenter-topenwide)); i<(int)(topencenter+topenwide); i++) {

      // how many this second?
      float fthis = MAX (0., topenmag * (0.5*cos( M_PI*((float)i-topencenter)/topenwide ) + 0.5));
      int nthis = (int)fthis;
      if ( fthis-(float)nthis > rand()/(float)RAND_MAX ) nthis++;

      // for each one
      float thisc[3];
      int globnumtries = 0;
      float frac_masked = 0.;
      for (int j=0; j<nthis; j++) {

         // timing
         block[icnt].tstart = (float)i + rand()/(float)RAND_MAX;
         block[icnt].tlen = 0.7 + 0.5*rand()/(float)RAND_MAX;
         block[icnt].mask = FALSE;

         // color
         (void) get_random_color (c,cnx,cny,thisc);
         block[icnt].r = thisc[0];
         block[icnt].g = thisc[1];
         block[icnt].b = thisc[2];

         // how much of the field is masked over?
         int num_masked = 0;
         for (int ix=0; ix<nx; ix++) for (int iy=0; iy<ny; iy++) if (masked[ix][iy]) num_masked++;
         frac_masked = (float)num_masked / (float)(nx*ny);

         // iterate until we uncover unique space
         int keepgoing = TRUE;
         int numtries = 0;
         while (keepgoing) {

            // location and size
            const int sizex = (int)sizebase + (int)(sizerand*rand()/(float)RAND_MAX);
            const int sizey = (int)sizebase + (int)(sizerand*rand()/(float)RAND_MAX);
            const int centerx = (float)nx * rand()/(float)RAND_MAX;
            const int centery = (float)ny * rand()/(float)RAND_MAX;
            block[icnt].startx = MAX(0,centerx-sizex/2);
            block[icnt].endx = MIN(nx,centerx+sizex/2);
            block[icnt].starty = MAX(0,centery-sizey/2);
            block[icnt].endy = MIN(ny,centery+sizey/2);

            // grow by 1 pixel if this leaves any 1-pixel mask walls

            // check bottom row
            int keeptrying = TRUE;
            while (keeptrying && (block[icnt].starty > 1)) {
               int isbad = FALSE;
               for (int ix=block[icnt].startx; ix<block[icnt].endx; ix++) {
                  if (masked[ix][block[icnt].starty-1] == TRUE &&
                      masked[ix][block[icnt].starty-2] == FALSE) {
                     isbad = TRUE;
                  }
               }
               if (isbad) {
                  block[icnt].starty--;
               } else {
                  keeptrying = FALSE;
               }
            }

            // check top row
            keeptrying = TRUE;
            while (keeptrying && (block[icnt].endy < ny-1)) {
               int isbad = FALSE;
               for (int ix=block[icnt].startx; ix<block[icnt].endx; ix++) {
                  if (masked[ix][block[icnt].endy] == TRUE &&
                      masked[ix][block[icnt].endy+1] == FALSE) {
                     isbad = TRUE;
                  }
               }
               if (isbad) {
                  block[icnt].endy++;
               } else {
                  keeptrying = FALSE;
               }
            }

            // check leftmost side
            keeptrying = TRUE;
            while (keeptrying && (block[icnt].startx > 1)) {
               int isbad = FALSE;
               for (int iy=block[icnt].starty; iy<block[icnt].endy; iy++) {
                  if (masked[block[icnt].startx-1][iy] == TRUE &&
                      masked[block[icnt].startx-2][iy] == FALSE) {
                     isbad = TRUE;
                  }
               }
               if (isbad) {
                  block[icnt].startx--;
               } else {
                  keeptrying = FALSE;
               }
            }

            // check rightmost side
            keeptrying = TRUE;
            while (keeptrying && (block[icnt].endx < nx-1)) {
               int isbad = FALSE;
               for (int iy=block[icnt].starty; iy<block[icnt].endy; iy++) {
                  if (masked[block[icnt].endx][iy] == TRUE &&
                      masked[block[icnt].endx+1][iy] == FALSE) {
                     isbad = TRUE;
                  }
               }
               if (isbad) {
                  block[icnt].endx++;
               } else {
                  keeptrying = FALSE;
               }
            }

            // test to see if it opens up enough
            int num_opened = 0;
            int num_pixels = 0;
            for (int ix=block[icnt].startx; ix<block[icnt].endx; ix++) {
               for (int iy=block[icnt].starty; iy<block[icnt].endy; iy++) {
                  num_pixels++;
                  if (masked[ix][iy]) num_opened++;
               }
            }
            float frac_opened = (float)num_opened / (float)num_pixels;

            numtries++;

            // when we get near the end, make sure we always open masked pixels
            if (frac_opened > frac_masked*frac_masked) keepgoing = FALSE;

            // but early on, try to adjoin existing open areas
            if (num_opened == num_pixels && rand()/(float)RAND_MAX < (float)icnt/30.) keepgoing = TRUE;

            // and don't try too hard
            if (numtries > 1000) keepgoing = FALSE;
         }
         globnumtries += numtries;

         // unmask those pixels
         for (int ix=block[icnt].startx; ix<block[icnt].endx; ix++) {
            for (int iy=block[icnt].starty; iy<block[icnt].endy; iy++) {
               masked[ix][iy] = FALSE;
            }
         }

         icnt++;
         if (icnt > MAXBLOCKS) break;
      }
      fprintf(stderr,"adding %d blocks at time %d, required %d tries, masked frac is %g\n",nthis,i,globnumtries,frac_masked); fflush(stderr);
   }
   const int iopened = icnt;

   // blocks that close off the space (completely)
   for (int i=(int)(tmaskcenter-tmaskwide); i<MIN(stop_time,(int)(tmaskcenter+tmaskwide)); i++) {

      // how many this second?
      float fthis = MAX (0., tmaskmag * (0.5*cos( M_PI*((float)i-tmaskcenter)/tmaskwide ) + 0.5));
      int nthis = (int)fthis;
      if ( fthis-(float)nthis > rand()/(float)RAND_MAX ) nthis++;

      // for each one
      float frac_masked = 0.;
      int globnumtries = 0;
      for (int j=0; j<nthis; j++) {

         // timing
         block[icnt].tstart = (float)i + rand()/(float)RAND_MAX;
         block[icnt].tlen = 0.5 + 0.4*rand()/(float)RAND_MAX;
         block[icnt].mask = TRUE;

         // color
         block[icnt].r = 0.;
         block[icnt].g = 0.;
         block[icnt].b = 0.;

         // how much of the field is masked over?
         int num_masked = 0;
         for (int ix=0; ix<nx; ix++) for (int iy=0; iy<ny; iy++) if (masked[ix][iy]) num_masked++;
         frac_masked = (float)num_masked / (float)(nx*ny);

         // iterate until we *cover* most unique unmasked space
         int keepgoing = TRUE;
         int numtries = 0;
         while (keepgoing) {

            // location and size
            const int sizex = (int)sizebase + (int)(sizerand*rand()/(float)RAND_MAX);
            const int sizey = (int)sizebase + (int)(sizerand*rand()/(float)RAND_MAX);
            const int centerx = (float)nx * rand()/(float)RAND_MAX;
            const int centery = (float)ny * rand()/(float)RAND_MAX;
            block[icnt].startx = MAX(0,centerx-sizex/2);
            block[icnt].endx = MIN(nx,centerx+sizex/2);
            block[icnt].starty = MAX(0,centery-sizey/2);
            block[icnt].endy = MIN(ny,centery+sizey/2);

            // test to see if it closes off enough
            int num_closed = 0;
            int num_pixels = 0;
            for (int ix=block[icnt].startx; ix<block[icnt].endx; ix++) {
               for (int iy=block[icnt].starty; iy<block[icnt].endy; iy++) {
                  num_pixels++;
                  if (!masked[ix][iy]) num_closed++;
               }
            }
            float frac_closed = (float)num_closed / (float)num_pixels;

            numtries++;

            // make sure we continue to mask off open areas
            if (frac_closed > (1.-frac_masked)*(1.-frac_masked)) keepgoing = FALSE;

            // but try to adjoin existing masked areas
            if (num_closed == num_pixels && rand()/(float)RAND_MAX < (float)(icnt-iopened)/30.) keepgoing = TRUE;

            // and don't try too hard
            if (numtries > 1000) keepgoing = FALSE;
         }
         globnumtries += numtries;

         // remask those pixels
         for (int ix=block[icnt].startx; ix<block[icnt].endx; ix++) {
            for (int iy=block[icnt].starty; iy<block[icnt].endy; iy++) {
               masked[ix][iy] = TRUE;
            }
         }

         icnt++;
         if (icnt > MAXBLOCKS) break;
      }
      fprintf(stderr,"masking %d blocks at time %d, required %d tries, masked frac is %g\n",nthis,i,globnumtries,frac_masked); fflush(stderr);
   }


   nblocks = icnt;

   free_2d_array_i(masked);
   free_2d_array_f(c[0]);
   free_2d_array_f(c[1]);
   free_2d_array_f(c[2]);

   return;
}

void update_mask_with_blocks_1 (float **mask, float **a[MAX_SCALARS],
                                int nx, int ny, float simtime) {

   // first test: just close off one random rectangle every frame
   const int sizex = 20 + (int)(100.0*rand()/(float)RAND_MAX);
   const int sizey = 20 + (int)(100.0*rand()/(float)RAND_MAX);
   const int centerx = (float)nx * rand()/(float)RAND_MAX;
   const int centery = (float)ny * rand()/(float)RAND_MAX;
   const int startx = MAX(0,centerx-sizex/2);
   const int endx = MIN(nx,centerx+sizex/2);
   const int starty = MAX(0,centery-sizey/2);
   const int endy = MIN(ny,centery+sizey/2);

   // close off the mask
   for (int ix=startx; ix<endx; ix++) {
      for (int iy=starty; iy<endy; iy++) {
         // use 0.0 for solid, 1.0 for open
         //mask[ix][iy] = 0.0;
      }
   }

   // and give it a random color
   const float r = rand()/(float)RAND_MAX;
   const float g = rand()/(float)RAND_MAX;
   const float b = rand()/(float)RAND_MAX;
   for (int ix=startx; ix<endx; ix++) {
      for (int iy=starty; iy<endy; iy++) {
         // use 0.0 for solid, 1.0 for open
         a[RR][ix][iy] = r;
         a[GG][ix][iy] = g;
         a[BB][ix][iy] = b;
         a[SF][ix][iy] = 0.3*r + 0.7*g + 0.1*b;
      }
   }

   return;
}

void update_mask_with_blocks_2 (float **mask, float **a[MAX_SCALARS],
                                int nx, int ny, float simtime, float dt) {

   // close off the mask
   //for (int ix=startx; ix<endx; ix++) {
   //   for (int iy=starty; iy<endy; iy++) {
         // use 0.0 for solid, 1.0 for open
         //mask[ix][iy] = 0.0;
   //   }
   //}

   // loop over all blocks
   int nthis = 0;
   for (int i=0; i<nblocks; i++) {
      // is this block fading in now?
      const float trel = (simtime - block[i].tstart) / block[i].tlen;
      if (trel > 0. && trel < 1.) {

         // weight fraction
         //float cratio = (1. - trel) / (1. - trel + dt/block[i].tlen);
         float cratio = (1. - trel - dt/block[i].tlen) / (1. - trel);
         if (cratio < 0.) cratio = 0.;
         //if (simtime+dt > block[i].tstart+block[i].tlen) cratio = 0.;
         fprintf(stderr,"cratio is %g\n",cratio); fflush(stderr);

         // color in the block
         for (int ix=block[i].startx; ix<block[i].endx; ix++) {
            for (int iy=block[i].starty; iy<block[i].endy; iy++) {
               // use 0.0 for solid, 1.0 for open
               a[RR][ix][iy] = cratio*a[RR][ix][iy] + (1.-cratio)*block[i].r;
               a[GG][ix][iy] = cratio*a[GG][ix][iy] + (1.-cratio)*block[i].g;
               a[BB][ix][iy] = cratio*a[BB][ix][iy] + (1.-cratio)*block[i].b;
               a[SF][ix][iy] = 0.3*a[RR][ix][iy] + 0.7*a[GG][ix][iy] + 0.1*a[BB][ix][iy];
               // flip immediately to new mask
               //if (block[i].mask) mask[ix][iy] = 0.;
               //else mask[ix][iy] = 1.;
               // smoothly shift to new mask
               if (block[i].mask) mask[ix][iy] = cratio*mask[ix][iy];
               else mask[ix][iy] = cratio*mask[ix][iy] + (1.-cratio);
            }
         }
         nthis++;
      }
   }
   fprintf(stderr,"drawing %d blocks this step time\n",nthis); fflush(stderr);

   return;
}

//
// Read, store, and overlay a second mask over the existing mask
//
void overlay_mask (int nx, int ny, float** mask) {

   static int first_time = TRUE;
   static float** overlay = NULL;
   //char maskfilename[255] = "burn_it_down_edges5_1025.png";
   char maskfilename[255] = "burn_it_down_edges3_2049.png";

   // if we're not using a mask, just ignore this
   if (mask == NULL) return;

   // read in the mask file
   if (first_time) {
      overlay = allocate_2d_array_f(nx,ny);
      // read grayscale PNG of exactly nx by ny resolution
      // 0.0 = black = masked, 1.0 = white = open
      read_png(maskfilename,nx,ny,FALSE,FALSE,1.0,FALSE,
         overlay,0.0,1.0,NULL,0.0,1.0,NULL,0.0,1.0);
      // do not normalize mask
      // never re-load or reallocate
      first_time = FALSE;
   }

   // lay it over the existing mask
   for (int ix=0; ix<nx; ix++) for (int iy=0; iy<ny; iy++)
      mask[ix][iy] = MIN(overlay[ix][iy], mask[ix][iy]);

   return;
}
