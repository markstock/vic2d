//
// Function to be passed into Fortran and called there to set wall BCs
//
// dp/dx + alfxa*p(xa,y) = gbdxa(y)
//
void funcbndyc (int* kbdy, float* xory, float* alfa, float* gbdy)
{

   static float vl = 0.0;
   static float vr = 0.0;
   static float ub = 0.0;
   static float ut = 0.0;

   // if kbdy is negative, we are calling this routine to save a boundary velocity
   if (*kbdy < 0 && *kbdy > -5) {
      if (*kbdy == -3) {
         // left (xs) wall
         vl = (*gbdy);
      } else if (*kbdy == -4) {
         // right (xf) wall
         vr = (*gbdy);
      } else if (*kbdy == -1) {
         // bottom (ys) wall
         ub = (*gbdy);
      } else {
         // top (yf) wall
         ut = (*gbdy);
      }
   }

   // for driven cavity sims, alfa is always 0
   (*alfa) = 0.0;
   (*gbdy) = 0.0;

   if (*kbdy == 3) {
      // left (xs) wall
      (*gbdy) = vl;
   } else if (*kbdy == 4) {
      // right (xf) wall
      (*gbdy) = vr;
   } else if (*kbdy == 1) {
      // bottom (ys) wall
      (*gbdy) = ub;
   } else {
      // top (yf) wall
      (*gbdy) = ut;
   }

   return;
}

// call cfx(x,cxx,cx,cex)
void funccfx(float* x, float* cxx, float* cx, float* cex) {
   (*cxx) = 1.0;
   (*cx)  = 0.0;
   (*cex) = 0.0;
   return;
}

void funccfy(float* y, float* cyy, float* cy, float* cey) {
   (*cyy) = 1.0;
   (*cy)  = 0.0;
   (*cey) = 0.0;
   return;
}

