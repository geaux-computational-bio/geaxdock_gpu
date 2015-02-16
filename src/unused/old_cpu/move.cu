#include <cmath>
#include <cstdio>

#include "dock.h"


// move the ligand around the pocket center
// calculate the center coorninate

void
Move (Ligand *lig, float tscale, float rscale)
{
  float tra_x = tscale * MyRand ();
  float tra_y = tscale * MyRand ();
  float tra_z = tscale * MyRand ();

  // euler angles
  float rot_x = rscale * MyRand ();
  float rot_y = rscale * MyRand ();
  float rot_z = rscale * MyRand ();

  // http://en.wikipedia.org/wiki/Euler_angles
  // http://upload.wikimedia.org/math/e/9/c/e9cf817bce9c1780216921cd93233459.png
  // http://upload.wikimedia.org/math/f/4/e/f4e55dc2c9581007648967d29b15121e.png
  // r = rz ry rx
  float s1 = sinf (rot_x);
  float c1 = cosf (rot_x);
  float s2 = sinf (rot_y);
  float c2 = cosf (rot_y);
  float s3 = sinf (rot_z);
  float c3 = cosf (rot_z);

  float r[3][3];
  r[0][0] = c1 * c2;
  r[0][1] = c1 * s2 * s3 - c3 * s1;
  r[0][2] = s1 * s3 + c1 * c3 * s2;
  r[1][0] = c2 * s1;
  r[1][1] = c1 * c3 + s1 * s2 * s3;
  r[1][2] = c3 * s1 * s2 - c1 * s3;
  r[2][0] = -1 * s2;
  r[2][1] = c2 * s3;
  r[2][2] = c2 * c3;



  int lig_idx = 0;
  Ligand *mylig = &lig[lig_idx];
  int n = mylig->lna;		// ligand point total number
  int track = mylig->track;
  Lig_Coord *coord_new = &mylig->coord[track];
  Lig_Coord *coord_old = &mylig->coord[!track];

  // center
  float cx = 0.0f;
  float cy = 0.0f;
  float cz = 0.0f;

  // iterate through all ligand atoms
  for (int atom = 0; atom < n; atom++) {
    float x = coord_old->x[atom];
    float y = coord_old->y[atom];
    float z = coord_old->z[atom];

    // rotation and translation
    float newx = r[0][0] * x + r[0][1] * y + r[0][2] * z + tra_x;
    float newy = r[1][0] * x + r[1][1] * y + r[1][2] * z + tra_y;
    float newz = r[2][0] * x + r[2][1] * y + r[2][2] * z + tra_z;

    // save to coord_new
    coord_new->x[atom] = newx;
    coord_new->y[atom] = newy;
    coord_new->z[atom] = newz;

    // accumulate center
    cx += newx;
    cy += newy;
    cz += newz;
  }

  // save center coornidate
  coord_new->center[0] = cx / n;
  coord_new->center[1] = cy / n;
  coord_new->center[2] = cz / n;
}

