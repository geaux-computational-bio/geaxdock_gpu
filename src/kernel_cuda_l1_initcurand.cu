// initialize curand status
// CURAND_Library.pdf, pp21

__global__ void
InitCurand_d ()
{
  const int gidx =
    blockDim.x * blockDim.y * blockIdx.x +
    blockDim.x * threadIdx.y +
    threadIdx.x;

  curand_init (seed_dc, gidx, 0, &curandstate_dc[gidx]);

  // seed, subsequence, offset, gpuseed
  // skipahead(100000, &curandstate_dc[gidx]);
}

