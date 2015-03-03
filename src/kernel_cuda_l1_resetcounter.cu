/*
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include <cuda.h>

#include "dock.h"
#include "gpu.cuh"
*/

__global__ void
ResetCounter_d (const int rep_begin, const int rep_end)
{
  const int bidx = blockDim.x * threadIdx.y + threadIdx.x;      // within a TB

  for (int offset = rep_begin; offset <= rep_end; offset += GD) {
    const int myreplica = offset + blockIdx.x;
    if (myreplica <= rep_end)
      if (bidx == 0)
	ligrecord_dc[myreplica - rep_begin].next_ptr = 0;

  }

}

 
