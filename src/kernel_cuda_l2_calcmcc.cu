/*
  #include <cmath>
  #include <cstdio>

  #include <cuda.h>

  #include "dock.h"
  #include "gpu.cuh"
*/

/*

__device__ void
serial_InitRefMatrix_d (const int bidx, Ligand * __restrict__ mylig, const Protein * __restrict__ myprt)
{


  if (bidx == 0) {


  // lig loop, ~30
    for (int l = 0; l < lna_dc; ++l) {
      const int lig_t = mylig->t[l];


      // prt loop, ~300
      for (int p = 0; p < pnp_dc; ++p) {
	  const int prt_t = myprt->t[p];

	  const float dx = mylig->coord_new.x[l] - myprt->x[p];
	  const float dy = mylig->coord_new.y[l] - myprt->y[p];
	  const float dz = mylig->coord_new.z[l] - myprt->z[p];
	  const float dst = sqrtf (dx * dx + dy * dy + dz * dz);

	  const float pmf0 = enepara_dc->pmf0[lig_t][prt_t];
	  ref_matrix_dc->matrix[l][p] = (dst <= pmf0);
	  //ref_matrix_dc->matrix[l][p] = dst;


	  //printf("init %f %d %d\n", ref_matrix_dc->matrix[l][p], l, p);
	  
      }			// prt loop
  }				// lig loop

    // printf ("init lna %d pnp %d pnk %d pos %d\n", lna_dc, pnp_dc, pnk_dc, pos_dc);

  }

}

*/





__device__ void
InitRefMatrix_d (const int bidx, Ligand * __restrict__ mylig, const Protein * __restrict__ myprt)
{
  // lig loop, ~30
  for (int i = 0; i < lna_dc; i += blockDim.y) {
    const int l = i + threadIdx.y;
    if (l < lna_dc) {
      const int lig_t = mylig->t[l];

      // prt loop, ~300
      for (int j = 0; j < pnp_dc; j += blockDim.x) {
	const int p = j + threadIdx.x;
	if (p < pnp_dc) {
	  const int prt_t = myprt->t[p];

	  const float dx = mylig->coord_new.x[l] - myprt->x[p];
	  const float dy = mylig->coord_new.y[l] - myprt->y[p];
	  const float dz = mylig->coord_new.z[l] - myprt->z[p];
	  const float dst = sqrtf (dx * dx + dy * dy + dz * dz);

	  const float pmf0 = enepara_dc->pmf0[lig_t][prt_t];
	  ref_matrix_dc->matrix[l][p] = (dst <= pmf0);
	}
      }				// prt loop
    }
  }				// lig loop
}



/*

__device__ void
serial_CalcMcc_d (const int bidx, Ligand * __restrict__ mylig, const Protein * __restrict__ myprt)
{

  if (bidx == 0) {
    float tp0 = 0.0f;
    float tn0 = 0.0f;
    float fp0 = 0.0f;
    float fn0 = 0.0f;

    // lig loop, ~30
    for (int l = 0; l < lna_dc; ++l) {
      const int lig_t = mylig->t[l];

      // prt loop, ~300
      for (int p = 0; p < pnp_dc; ++p) {

	//printf ("%d %d, ", l, p);

	const int prt_t = myprt->t[p];

	const float dx = mylig->coord_new.x[l] - myprt->x[p];
	const float dy = mylig->coord_new.y[l] - myprt->y[p];
	const float dz = mylig->coord_new.z[l] - myprt->z[p];
	const float dst = sqrtf (dx * dx + dy * dy + dz * dz);

	const float pmf0 = enepara_dc->pmf0[lig_t][prt_t];
	const int ref_val = ref_matrix_dc->matrix[l][p];

	//if (dst != ref_val)
	//printf("calc %f %f %d %d\n", dst, ref_val, l, p);

	tp0 += (ref_val == 1 && dst <= pmf0);
	fn0 += (ref_val == 1 && dst > pmf0);
	fp0 += (ref_val == 0 && dst <= pmf0);
	tn0 += (ref_val == 0 && dst > pmf0);

      }				// prt loop
    }				// lig loop

    const float dividend = sqrtf ((tp0 + fp0) * (tp0 + fn0) * (tn0 + fp0) * (tn0 + fn0));

    if (dividend != 0)
      mylig->energy_new.cms = (tp0 * tn0 - fp0 * fn0) / dividend;
    else
      mylig->energy_new.cms = CMCC_INVALID_VAL;

    // printf("%f %f %f %f %f\n", tp0, fn0, fp0, tn0, dividend);
    // printf("cms %f\n", mylig->energy_new.cms);

    //printf ("calc lna %d pnp %d pnk %d pos %d\n", lna_dc, pnp_dc, pnk_dc, pos_dc);

  }

}
*/



__device__ void
CalcMcc_d (const int bidx, Ligand * __restrict__ mylig, const Protein * __restrict__ myprt)
{
  // reduce
  __shared__ int tp[TperB];
  __shared__ int fn[TperB];
  __shared__ int fp[TperB];
  __shared__ int tn[TperB];
  tp[bidx] = 0;
  fn[bidx] = 0;
  fp[bidx] = 0;
  tn[bidx] = 0;
  __syncthreads ();

  // lig loop, ~30
  for (int i = 0; i < lna_dc; i += blockDim.y) {
    const int l = i + threadIdx.y;
    if (l < lna_dc) {
      const int lig_t = mylig->t[l];

      // prt loop, ~300
      for (int j = 0; j < pnp_dc; j += blockDim.x) {
	const int p = j + threadIdx.x;
	if (p < pnp_dc) {

	  const int prt_t = myprt->t[p];

	  const float dx = mylig->coord_new.x[l] - myprt->x[p];
	  const float dy = mylig->coord_new.y[l] - myprt->y[p];
	  const float dz = mylig->coord_new.z[l] - myprt->z[p];
	  const float dst = sqrtf (dx * dx + dy * dy + dz * dz);

	  const float pmf0 = enepara_dc->pmf0[lig_t][prt_t];
	  const int ref_val = ref_matrix_dc->matrix[l][p];

	  tp[bidx] += (ref_val == 1 && dst <= pmf0);
	  fn[bidx] += (ref_val == 1 && dst > pmf0);
	  fp[bidx] += (ref_val == 0 && dst <= pmf0);
	  tn[bidx] += (ref_val == 0 && dst > pmf0);
	}
      }				// prt loop
    }
  }				// lig loop

  SumReduction_int_1D_4_d (bidx, tp, fn, fp, tn);

  if (bidx == 0) {
    const float tp0 = (float) tp[0];
    const float fn0 = (float) fn[0];
    const float fp0 = (float) fp[0];
    const float tn0 = (float) tn[0];


    const float dividend = sqrtf ((tp0 + fp0) * (tp0 + fn0) * (tn0 + fp0) * (tn0 + fn0));

    if (dividend != 0)
      mylig->energy_new.cms = (tp0 * tn0 - fp0 * fn0) / dividend;
    else
      mylig->energy_new.cms = CMCC_INVALID_VAL;

    // printf("%f %f %f %f %f\n", tp0, fn0, fp0, tn0, dividend);
    // printf("cms %f\n", mylig->energy_new.cms);

  }

}
