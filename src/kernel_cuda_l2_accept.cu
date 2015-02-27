/*
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include <cuda.h>

#include "dock.h"
#include "gpu.cuh"
*/

__forceinline__
__device__ void
Accept_d (const int bidx, Ligand * __restrict__ mylig, const float mybeta, const int myreplica)
{
  __shared__ int is_accept;

  if (bidx == 0) {    
#if IS_FORCE_TO_ACCEPT == 1
      is_accept = 1;
#elif IS_FORCE_TO_ACCEPT == 0
      const float delta_energy = mylig->energy_new.e[MAXWEI - 1] - mylig->energy_old.e[MAXWEI -1];
      is_accept = (MyRand_d () < expf (delta_energy * mybeta));  // mybeta is less than zero
      // printf ("delta_energy: %.8f\n", delta_energy);
      // printf("is_accept: %d\n", is_accept);
      // printf("Myrand_d: %f\n", MyRand_d());
      // printf("prob: %.32f\n", expf (delta_energy * mybeta));
      // printf("mybeta: %.20f\n", mybeta);
#endif
    //printf ("beta[%d] = %f\n", myreplica, 1.0f / mybeta);
    mylig->is_move_accepted = is_accept;
  }


  __syncthreads ();

  if (is_accept == 1) {
    if (bidx < 6)
      mylig->movematrix_old[bidx] = mylig->movematrix_new[bidx];
    if (bidx == 0) {
	  mylig->energy_old = mylig->energy_new;
    }

    /*
    if (bidx == 0 && myreplica == 0 && is_accept == 1)
      printf ("accept_d: accepted %d\n", acs_mc_dc[myreplica]);
    */

  }


}

 
