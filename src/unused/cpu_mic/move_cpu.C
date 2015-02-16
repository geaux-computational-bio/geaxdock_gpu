/*
#include <cmath>
#include <cstdio>

#include "dock.h"
#include "gpu.cuh"
*/
#include "debug.h"


// move the ligand around the pocket center

void
Move_d (Ligand * __restrict__ mylig)
{
  float movematrix_new[6]; // translation x y z, rotation x y z
  float rot[3][3]; // rotz roty rotx


  for (int bidx = 0; bidx < 6; ++bidx) {
    // range of perturbation should be [-1 , 1]
    float perturbation = 2.0f * MyRand_d() - 1.0f;  
    movematrix_new[bidx] = move_scale_dc[bidx] * perturbation + mylig->movematrix_old[bidx];
    mylig->movematrix_new[bidx] = movematrix_new[bidx];
  }


  // http://en.wikipedia.org/wiki/Euler_angles
  // http://upload.wikimedia.org/math/e/9/c/e9cf817bce9c1780216921cd93233459.png
  // http://upload.wikimedia.org/math/f/4/e/f4e55dc2c9581007648967d29b15121e.png
  const float s1 = sinf (movematrix_new[3]);
  const float c1 = cosf (movematrix_new[3]);
  const float s2 = sinf (movematrix_new[4]);
  const float c2 = cosf (movematrix_new[4]);
  const float s3 = sinf (movematrix_new[5]);
  const float c3 = cosf (movematrix_new[5]);
  
  rot[0][0] = c1 * c2;
  rot[0][1] = c1 * s2 * s3 - c3 * s1;
  rot[0][2] = s1 * s3 + c1 * c3 * s2;
  rot[1][0] = c2 * s1;
  rot[1][1] = c1 * c3 + s1 * s2 * s3;
  rot[1][2] = c3 * s1 * s2 - c1 * s3;
  rot[2][0] = -1 * s2;
  rot[2][1] = c2 * s3;
  rot[2][2] = c2 * c3;



  LigCoord *coord_new = &mylig->coord_new;
  LigCoord *coord_orig = &mylig->coord_orig;


  const float cx = coord_orig->center[0];
  const float cy = coord_orig->center[1];
  const float cz = coord_orig->center[2];
  
  // iterate through all ligand residues
  // rotation and translation, and apply coordinate system transformation
  for (int l = 0; l < lna_dc; ++l) {
    float x = coord_orig->x[l];
    float y = coord_orig->y[l];
    float z = coord_orig->z[l];
    coord_new->x[l] = rot[0][0] * x + rot[0][1] * y + rot[0][2] * z + movematrix_new[0] + cx;
    coord_new->y[l] = rot[1][0] * x + rot[1][1] * y + rot[1][2] * z + movematrix_new[1] + cy;
    coord_new->z[l] = rot[2][0] * x + rot[2][1] * y + rot[2][2] * z + movematrix_new[2] + cz;
  }
  
  
  for (int i = 0; i < 3; ++i) { 
    coord_new->center[i] = coord_orig->center[i] + movematrix_new[i];
  }
  
}


