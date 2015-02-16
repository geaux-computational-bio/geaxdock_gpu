__device__ void
CalcRmsd_d (const int bidx, Ligand * __restrict__ mylig)
{
  
  __shared__ float distance_square[BDx];

  if (threadIdx.y == 0) {
    distance_square[threadIdx.x] = 0.0f;
    LigCoord *coord_new = &mylig->coord_new;
    LigCoord *coord_orig = &mylig->coord_orig;
    
    const float orig_cx = coord_orig->center[0];
    const float orig_cy = coord_orig->center[1];
    const float orig_cz = coord_orig->center[2];

    for (int l = threadIdx.x; l < lna_dc; l += BDx) {
      float d_x = coord_new->x[l] - (coord_orig->x[l] + orig_cx);
      float d_y = coord_new->y[l] - (coord_orig->y[l] + orig_cy);
      float d_z = coord_new->z[l] - (coord_orig->z[l] + orig_cz);

      distance_square[l] += d_x * d_x + d_y * d_y + d_z * d_z;
    }
  }


  // reduce for the sum of distance squares
  __syncthreads (); 
  for (int stride = BDx / 2; stride >= 1; stride >>= 1) {
    if (bidx < stride)
      distance_square[bidx] += distance_square[bidx + stride];
    __syncthreads (); 
  }

  if (bidx == 0)
    mylig->energy_new.rmsd = sqrtf (distance_square[0] / lna_dc);
}
