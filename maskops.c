#include "vicmoc.h"

typedef struct node_record *node_ptr;
typedef struct node_record {
   int index;
   node_ptr next;
} NODE;

void update_mask_with_blocks (float **mask, float **a[MAX_SCALARS],
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
}
