/*
#include <cstdlib>
#include <cstdio>

#include "dock.h"
#include "gpu.cuh"

#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
*/





__device__ void
InitAcs_d (const int bidx)
{
  if (blockIdx.x == 0) {
    for (int i = bidx; i < MAXREP; i += TperB) {
      acs_mc_dc[i] = 0;
      acs_temp_exchg_dc[i] = 0;
    }
  }
}

__device__ void
InitLigRecord_d (const int bidx, const int myreplica, const int rep_begin)
{
  for (int s2s3 = 0; s2s3 < steps_per_dump_dc; ++s2s3) {
    LigRecordSingleStep *myrecord =
      &ligrecord_dc[myreplica - rep_begin].step[s2s3];
    myrecord->replica.idx_rep = 0;
    myrecord->replica.idx_prt = 0;
    myrecord->replica.idx_tmp = 0;
    myrecord->replica.idx_lig = 0;

    for (int i = 0; i < MAXWEI; ++i)
      myrecord->energy.e[i] = 0.0f;

    for (int i = 0; i < 6; ++i)
      myrecord->movematrix[i] = 0.0f;

    myrecord->step = 0;
  }

}


/*
__forceinline__
__device__ void
BackupLigCoord_d (const int bidx, Ligand *mylig)
{

  const LigCoord *src = &mylig->coord_old;
  LigCoord *dst = &mylig->coord_bkup;

  for (int atom = bidx; atom < lna_dc; atom += TperB) {
    dst->x[atom] = src->x[atom];
    dst->y[atom] = src->y[atom];
    dst->z[atom] = src->z[atom];
  }
  if (bidx < 3)
    dst->center[bidx] = src->center[bidx];

}
*/



__device__ void
RecordLigand_d (const int bidx, const int s1, const int s2s3,
		const int myreplica, const int rep_begin,
		const Ligand * mylig, const int is_record)
{
  /*
     if (bidx == 0) // && myreplica == 0)
     printf ("rep %d, iter %d, rep_begin %d, n_rep %d, idx_rep %d\n",
     myreplica, s2, rep_begin, n_rep_dc, replica_dc[myreplica].idx_rep);

     if (myreplica == 0) {
     PrintEnergy2_d (bidx, mylig, myreplica, s1 + s2s3, 2);
     }
   */

  if (is_record == 1 && bidx == 0) {
    const int next_ptr = ligrecord_dc[myreplica - rep_begin].next_ptr;
    ligrecord_dc[myreplica - rep_begin].next_ptr = next_ptr + 1;
    LigRecordSingleStep *myrecord = &ligrecord_dc[myreplica - rep_begin].step[next_ptr];

    myrecord->replica = replica_dc[myreplica];
    myrecord->energy = mylig->energy_old;
    for (int i = 0; i < 6; ++i)
      myrecord->movematrix[i] = mylig->movematrix_old[i];
    myrecord->step = s1 + s2s3;
  }

}






__forceinline__ __device__ float
MyRand_d ()
{
  const int gidx =
    blockDim.x * blockDim.y * blockIdx.x +
    blockDim.x * threadIdx.y + threadIdx.x;
  curandState myseed = curandstate_dc[gidx];
  float randdd = curand_uniform (&myseed);
  curandstate_dc[gidx] = myseed;

  return randdd;
}




/*
__forceinline__
__device__ int
Mininal_int_d (const int a, const int b)
{
 return a < b ? a : b;
}
*/





__forceinline__ __device__ void
SumReduction1D_d (const int bidx, float *a)
{
  __syncthreads ();

  for (int stride = TperB / 2; stride >= 1; stride >>= 1) {
    if (bidx < stride)
      a[bidx] += a[stride + bidx];
    __syncthreads ();
  }
}

__forceinline__ __device__ void
SumReduction_int_1D_4_d (const int bidx, int *a, int *b, int *c, int *d)
{
  __syncthreads ();

  for (int stride = TperB / 2; stride >= 1; stride >>= 1) {
    if (bidx < stride) {
      a[bidx] += a[stride + bidx];
      b[bidx] += b[stride + bidx];
      c[bidx] += c[stride + bidx];
      d[bidx] += d[stride + bidx];
    }
    __syncthreads ();
  }
}

__forceinline__ __device__ void
SumReduction1D_5_d (const int bidx, float *a, float *b, float *c, float *d,
		    float *e)
{
  __syncthreads ();

  for (int stride = TperB / 2; stride >= 1; stride >>= 1) {
    if (bidx < stride) {
      a[bidx] += a[stride + bidx];
      b[bidx] += b[stride + bidx];
      c[bidx] += c[stride + bidx];
      d[bidx] += d[stride + bidx];
      e[bidx] += e[stride + bidx];
    }
    __syncthreads ();
  }
}


__forceinline__ __device__ void
SumReduction2D_d (float a[BDy][BDx])
{
  __syncthreads ();

  for (int stride = BDx / 2; stride >= 1; stride >>= 1) {
    if (threadIdx.x < stride) {
      a[threadIdx.y][threadIdx.x] += a[threadIdx.y][stride + threadIdx.x];
    }
    __syncthreads ();
  }
}


__forceinline__ __device__ void
SumReduction2D_2_d (float a[BDy][BDx], float b[BDy][BDx])
{
  __syncthreads ();

  for (int stride = BDx / 2; stride >= 1; stride >>= 1) {
    if (threadIdx.x < stride) {
      a[threadIdx.y][threadIdx.x] += a[threadIdx.y][stride + threadIdx.x];
      b[threadIdx.y][threadIdx.x] += b[threadIdx.y][stride + threadIdx.x];
    }
    __syncthreads ();
  }
}


__forceinline__ __device__ float
NormPdf (float x, float loc, float scale)
{

  float norm_para, prob, pdf_val;

  norm_para = 1 / (scale * sqrt (2 * PI));
  prob = exp (0.f - (x - loc) * (x - loc) / (2 * scale * scale));

  pdf_val = norm_para * prob;

  return pdf_val;
}

__forceinline__ __device__ float
CauchyPdf (float x, float loc, float scale)
{
  float norm_para, prob, pdf_val;

  norm_para = 1 / (PI * scale);
  prob = 1 / (1 + ((x - loc) / scale) * ((x - loc) / scale));

  pdf_val = norm_para * prob;

  return pdf_val;
}


__forceinline__ __device__ float
LogisticPdf (float x, float loc, float scale)
{
  float norm_para, e_power, prob, pdf_val;

  norm_para = 1 / scale;
  e_power = exp (-(x - loc) / scale);
  prob = e_power / powf (1 + e_power, 2.0);

  pdf_val = norm_para * prob;

  return pdf_val;
}

__forceinline__ __device__ float
WaldPdf (float x, float loc, float scale)
{
  float norm_para, prob, pdf_val;

  float normed_x = (x - loc) / scale;

  norm_para = 1 / (sqrt (2 * PI * powf (normed_x, 3.0)) * scale);
  prob = exp (-pow (normed_x - 1, 2) / (2 * normed_x));

  if (normed_x < 0)
    pdf_val = 0.00000001f;
  else
    pdf_val = norm_para * prob;

  return pdf_val;
}

__forceinline__ __device__ float
LaplacePdf (float x, float loc, float scale)
{
  float normed_x, pdf_val;

  normed_x = fabs (x - loc) / scale;

  pdf_val = (1 / (2 * scale)) * exp (-normed_x);

  return pdf_val;
}
