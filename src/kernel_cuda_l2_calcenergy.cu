
/*
#include <cmath>
#include <cstdio>

#include <cuda.h>

#include "dock.h"
#include "gpu.cuh"
*/



/*
#define expf(a) (a)
#define powf(a,b) (a+b)
#define logf(a) (a)
#define sqrtf(a) (a)
*/



__device__ void
CalcEnergy_d (const int bidx, Ligand * __restrict__ mylig, const Protein * __restrict__ myprt)
{
  const float sqrt_2_pi_m1 = -1.0f / sqrtf (2.0f * PI);


  // reduce all points on the X-Y plate
  __shared__ float evdw[TperB]; // e[0]
  __shared__ float eele[TperB]; // e[1]
  __shared__ float epmf[TperB]; // e[2]
  __shared__ float epsp[TperB]; // e[3]
  __shared__ float ehdb[TperB]; // e[4]

  // reduce through only x axis
  __shared__ float a_val[BDy][BDx]; // reused by hpc, kde, lhm
  __shared__ float a_sz[BDy][BDx];

  __shared__ float ehpc[BDy]; // e[5]
  __shared__ float ekde[BDy]; // e[6]
  __shared__ float elhm[BDy]; // e[7]


  __shared__ float enepara_p1a[MAXTP2][MAXTP1];
  __shared__ float enepara_p2a[MAXTP2][MAXTP1];
  __shared__ float enepara_pmf0[MAXTP2][MAXTP1];
  __shared__ float enepara_pmf1[MAXTP2][MAXTP1];
  __shared__ float enepara_hdb0[MAXTP2][MAXTP1];
  __shared__ float enepara_hdb1[MAXTP2][MAXTP1];
  __shared__ float enepara_hpl0[MAXTP2];
  __shared__ float enepara_hpl1[MAXTP2];
  __shared__ float enepara_hpl2[MAXTP2];



  evdw[bidx] = 0.0f;
  eele[bidx] = 0.0f;
  epmf[bidx] = 0.0f;
  epsp[bidx] = 0.0f;
  ehdb[bidx] = 0.0f;

  if (bidx < BDy) {
    ehpc[bidx] = 0.0f;
    ekde[bidx] = 0.0f;
    elhm[bidx] = 0.0f;
  }



  for (int i = 0; i < MAXTP2; i += blockDim.y) {
    const int l = i + threadIdx.y;
    if (l < MAXTP2) {

      for (int j = 0; j < MAXTP1; j += blockDim.x) {
	const int p = j + threadIdx.x;
	if (p < MAXTP1) {
	  enepara_p1a[l][p] = enepara_dc->p1a[l][p];
	  enepara_p2a[l][p] = enepara_dc->p2a[l][p];
	  enepara_pmf0[l][p] = enepara_dc->pmf0[l][p];
	  enepara_pmf1[l][p] = enepara_dc->pmf1[l][p];
	  enepara_hdb0[l][p] = enepara_dc->hdb0[l][p];
	  enepara_hdb1[l][p] = enepara_dc->hdb1[l][p];
	}
      }

    }
  }

  for (int i = 0; i < MAXTP2; i += TperB) {
    const int l = i + bidx;
    if (l < MAXTP2) {
      enepara_hpl0[l] = enepara_dc->hpl0[l];
      enepara_hpl1[l] = enepara_dc->hpl1[l];
      enepara_hpl2[l] = enepara_dc->hpl2[l];
    }
  }



  __syncthreads ();





#if 1
  // lig loop, ~30
  for (int i = 0; i < lna_dc; i += blockDim.y) {
    a_val[threadIdx.y][threadIdx.x] = 0.0f;
    const int l = i + threadIdx.y;
    if (l < lna_dc) {
      const int lig_t = mylig->t[l];

      // prt loop, ~300
      for (int j = 0; j < pnp_dc; j += blockDim.x) {
	const int p = j + threadIdx.x;
	if (p < pnp_dc) {
	  
	  const int prt_t = CUDA_LDG_D (myprt->t[p]);

	  const float dx = mylig->coord_new.x[l] - CUDA_LDG_D (myprt->x[p]);
	  const float dy = mylig->coord_new.y[l] - CUDA_LDG_D (myprt->y[p]);
	  const float dz = mylig->coord_new.z[l] - CUDA_LDG_D (myprt->z[p]);
	  const float dst_pow2 = dx * dx + dy * dy + dz * dz;
	  const float dst_pow4 = dst_pow2 * dst_pow2;
	  const float dst = sqrtf (dst_pow2);


#if 1
	  /* hydrophobic potential */
	  if (CUDA_LDG_D (myprt->c0_and_d12_or_c2[p]) == 1 && dst_pow2 <= 81.0f) {
	    a_val[threadIdx.y][threadIdx.x] += CUDA_LDG_D (myprt->hpp[p]) *
	      (1.0f - (3.5f / 81.0f * dst_pow2 -
		       4.5f / 81.0f / 81.0f * dst_pow4 +
		       2.5f / 81.0f / 81.0f / 81.0f * dst_pow4 * dst_pow2 -
		       0.5f / 81.0f / 81.0f / 81.0f / 81.0f * dst_pow4 * dst_pow4));
	  }
#endif
      

#if 1
	  // p1a[MAXTP2][MAXTP1]
	  // p2a[MAXTP2][MAXTP1]

	  /* L-J potential */
	  const float p1 = enepara_p1a[lig_t][prt_t] / (dst_pow4 * dst_pow4 * dst);
	  const float p2 = enepara_p2a[lig_t][prt_t] / (dst_pow4 * dst_pow2);
	  const float p4 = p1 * enepara_lj0_dc * (1.0f + enepara_lj1_dc * dst_pow2) + 1.0f;
	  evdw[bidx] += (p1 - p2) / p4;
#endif




#if 1
	  /* electrostatic potential */
	  const float s1 = enepara_el1_dc * dst;
	  float g1;
	  if (s1 < 1)
	    g1 = enepara_el0_dc + enepara_a1_dc * s1 * s1 + enepara_b1_dc * s1 * s1 * s1;
	  else
	    g1 = 1.0f / s1;
	  eele[bidx] += CUDA_LDG_D (mylig->c[l]) * CUDA_LDG_D (myprt->ele[p]) * g1;
#endif

      
#if 1
	  // pmf0[MAXTP2][MAXTP1]
	  // pmf1[MAXTP2][MAXTP1]
	  // psp[MAXTP2][MAXPRO]

	  /* contact potential */
	  const float dst_minus_pmf0 = dst - enepara_pmf0[lig_t][prt_t];

	  epmf[bidx] +=
	    enepara_pmf1[lig_t][prt_t] /
	    (1.0f + expf ((-0.5f * dst + 6.0f) * dst_minus_pmf0));


	  /* pocket-specific potential */
	  // the senmatics do not match with the original program:
	  // if (found psp[][])
	  //   accumulate to epsp;
	  // else
	  //   do nothing
	  if (CUDA_LDG_D (myprt->c[p]) == 2 && dst_minus_pmf0 <= 0) {
	    const int i1 = CUDA_LDG_D (myprt->seq3r[p]);
	    epsp[bidx] += CUDA_LDG_D (psp_dc->psp[lig_t][i1]); // sparse matrix
	  }
#endif


#if 1
	  // hdb0[MAXTP2][MAXTP1]
	  // hdb1[MAXTP2][MAXTP1]



	  /* hydrogen bond potential */
	  const float hdb0 = enepara_hdb0[lig_t][prt_t];
	  if (hdb0 > 0.1f) {
	    const float hdb1 = enepara_hdb1[lig_t][prt_t];
	    const float hdb3 = (dst - hdb0) * hdb1;
	    ehdb[bidx] += sqrt_2_pi_m1 * hdb1 * expf (-0.5f * hdb3 * hdb3);
	  }
#endif

	} // if (p < pnp_dc)
      } // prt loop
    } // if (l < lna_dc)


#if 1
    // not performance critical
    // hpl0[MAXTP2]
    // hpl1[MAXTP2]
    // hpl2[MAXTP2]

    /* hydrophobic restraits*/
    SumReduction2D_d (a_val);
    // transpose may help improve the performance
    if (threadIdx.x == 0 && l < lna_dc) {
      const int lig_t = CUDA_LDG_D (mylig->t[l]);
      const float hpc2 = (a_val[threadIdx.y][0] - enepara_hpl0[lig_t]) / enepara_hpl1[lig_t];
      ehpc[threadIdx.y] += 0.5f * hpc2 * hpc2 - enepara_hpl2[lig_t];
    }
#endif



  } // lig loop
#endif


  SumReduction1D_5_d (bidx, evdw, eele, epmf, epsp, ehdb);


  if (bidx == 0) {
    float eehpc = 0.0f;
    for (int i = 0; i < BDy; ++i)
      eehpc += ehpc[i];
    ehpc[0] = eehpc;
  }




#if 1
  /* kde potential */

  // lig loop, ~30
  for (int i = 0; i < lna_dc; i += blockDim.y) {
    a_val[threadIdx.y][threadIdx.x] = 0.0f;
    a_sz[threadIdx.y][threadIdx.x] = 0.0f;
    const int l = i + threadIdx.y;
    if (l < lna_dc) {

      // kde loop, ~400
      for (int j = 0; j < pnk_dc; j += blockDim.x) {
	const int k = j + threadIdx.x;
	if (k < pnk_dc) {

	  if (CUDA_LDG_D (mylig->t[l]) == kde_dc->t[k]) {
	    const float dx = mylig->coord_new.x[l] - kde_dc->x[k];
	    const float dy = mylig->coord_new.y[l] - kde_dc->y[k];
	    const float dz = mylig->coord_new.z[l] - kde_dc->z[k];
	    const float kde_dst_pow2 = dx * dx + dy * dy + dz * dz;
	    a_val[threadIdx.y][threadIdx.x] += expf (enepara_kde2_dc * kde_dst_pow2);
	    a_sz[threadIdx.y][threadIdx.x] += 1.0f;
	  }

	} // if (k < pnk_dc)
      } // kde loop
    } // if (l < lna_dc)

    SumReduction2D_2_d (a_val, a_sz);

    if (threadIdx.x == 0 && l < lna_dc && a_sz[threadIdx.y][0] != 0.0f)
      ekde[threadIdx.y] += (a_val[threadIdx.y][0] / a_sz[threadIdx.y][0]);

  } // lig loop

  __syncthreads ();
  if (bidx == 0) {
    float eekde = 0.0f;
    for (int i = 0; i < BDy; ++i)
      eekde += ekde[i];
    eekde = eekde / enepara_kde3_dc;
    ekde[0] = eekde;
  }
  __syncthreads ();

#endif


#if 1

  /* position restraints */

  // lhm loop, ~11
  for (int i = 0; i < pos_dc; i += blockDim.y) {
    a_val[threadIdx.y][threadIdx.x] = 0.0f;
    a_sz[threadIdx.y][threadIdx.x] = 0.0f;
    const int m = i + threadIdx.y;

    if (m < pos_dc) {

    // lig loop, ~30
      for (int j = 0; j < lna_dc; j += blockDim.x) {
	const int l = j + threadIdx.x;
	if (l < lna_dc) {
	  const int lig_n = CUDA_LDG_D (mylig->n[l]) + 1;
	  if (mcs_dc[m].x[lig_n] != MCS_INVALID_COORD) {
	    const float dx = mylig->coord_new.x[l] - mcs_dc[m].x[lig_n];
	    const float dy = mylig->coord_new.y[l] - mcs_dc[m].y[lig_n];
	    const float dz = mylig->coord_new.z[l] - mcs_dc[m].z[lig_n];
	    a_val[threadIdx.y][threadIdx.x] += dx * dx + dy * dy + dz * dz;
	    a_sz[threadIdx.y][threadIdx.x] += 1.0f;
	  }
	} // if (l < lna_dc)
      } // lig loop

    } // if (m < pos_dc)

    SumReduction2D_2_d (a_val, a_sz);

    if (threadIdx.x == 0 && m < pos_dc) {
      elhm[threadIdx.y] +=
	mcs_dc[m].tcc *
	sqrtf (a_val[threadIdx.y][0] / a_sz[threadIdx.y][0]);
    }
  } // lhm loop

  __syncthreads ();
  if (bidx == 0) {
    float eelhm = 0.0f;
    for (int i = 0; i < BDy; ++i)
      eelhm += elhm[i];
    elhm[0] = eelhm;
  }
  __syncthreads ();

#endif





__shared__ float edst;
if (bidx == 0)
  edst = 9.9;
#if 0

  // energy edst e[8]
  __shared__ float edst;

  if (bidx == 0) {
    const float dx = mylig->coord_new.center[0] - myprt->pocket_center[0];
    const float dy = mylig->coord_new.center[1] - myprt->pocket_center[1];
    const float dz = mylig->coord_new.center[2] - myprt->pocket_center[2];
    dst = sqrtf (dx * dx + dy * dy + dz * dz);
  }
  __syncthreads ();

#endif


  Energy e;

  if (bidx == 0) {
    e.e[0] = evdw[0] / lna_dc; // 0 - vdw 
    e.e[1] = eele[0] / lna_dc; // 1 - ele
    e.e[2] = epmf[0] / lna_dc; // 2 - pmf (CP)
    e.e[3] = epsp[0] / lna_dc; // 3 - psp (PS CP)
    e.e[4] = ehdb[0] / lna_dc; // 4 - hdb (HB)
    e.e[5] = ehpc[0] / lna_dc; // 5 - hpc (HP)
    e.e[6] = ekde[0] / lna_dc; // 6 - kde (PHR)
    e.e[7] = logf (elhm[0] / pos_dc); // 7 - lhm (MCS)
    e.e[8] = edst; // 8 - dst (DST)

    // normalization
    for (int i = 0; i < MAXWEI - 1; ++i)
      e.e[i] = CUDA_LDG_D (enepara_dc->a_para[i]) * e.e[i] + CUDA_LDG_D (enepara_dc->b_para[i]);
  }


  // calculate the total energy from energy terms
  CombineEnergy_d (bidx, &e);


  if (bidx == 0)
    mylig->energy_new = e;
}


