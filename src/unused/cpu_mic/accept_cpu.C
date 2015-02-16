/*
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include <cuda.h>

#include "dock.h"
#include "gpu.cuh"
*/


void
Accept_d (Ligand * __restrict__ mylig, const float mybeta, const int myreplica)
{
  int is_accept;

  if (is_random_dc == 0) {
    is_accept = 1;
  }
  else {
    const float delta_energy = mylig->energy_new.e[MAXWEI - 1] - mylig->energy_old.e[MAXWEI -1];
    is_accept = (MyRand_d () < expf (delta_energy * mybeta));  // mybeta is less than zero
    // printf("is_accept: %d\n", is_accept);
    // printf("Myrand_d: %f\n", MyRand_d());
    // printf("prob: %f\n", expf (delta_energy * mybeta));
  }

  acs_dc[myreplica] += is_accept;


  if (is_accept == 1) {
    for (int i = 0; i < 6; ++i)
      mylig->movematrix_old[i] = mylig->movematrix_new[i];
    for (int i = 0; i < MAXWEI; ++i)
      mylig->energy_old.e[i] = mylig->energy_new.e[i];
  }
}

 
