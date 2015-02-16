
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




void
CalcEnergy_d (Ligand * __restrict__ mylig, const Protein * myprt)
{
  // reduce all points on the X-Y plate
  float evdw = 0.0f; // e[0]
  float eele = 0.0f; // e[1]
  float epmf = 0.0f; // e[2]
  float epsp = 0.0f; // e[3]
  float ehdb = 0.0f; // e[4]
  float ehpc = 0.0f; // e[5]
  float ekde = 0.0f; // e[6]
  float elhm = 0.0f; // e[7]



/* icc -mmic -vec-report3: loop was not vectorized: not inner loop */
  //#pragma omp parallel for reduction(+:evdw, eele, epmf, epsp, ehdb, ehpc)
  // lig loop, ~30
  for (int l = 0; l < lna_dc; ++l) {
    float hpc1 = 0.0f;
    const int lig_t = mylig->lig_point[l].t;

/* icc -mmic -vec-report3: LOOP WAS VECTORIZED */
    // prt loop, ~300
    for (int p = 0; p < pnp_dc; ++p) {

      const int prt_t = myprt->prt_point[p].t;

      const float dx = mylig->coord_new.x[l] - myprt->prt_point[p].x;
      const float dy = mylig->coord_new.y[l] - myprt->prt_point[p].y;
      const float dz = mylig->coord_new.z[l] - myprt->prt_point[p].z;
      const float dst_pow2 = dx * dx + dy * dy + dz * dz;
      const float dst_pow4 = dst_pow2 * dst_pow2;
      const float dst = sqrtf (dst_pow2);
      

      /* hydrophobic potential */
      if (myprt->prt_point[p].c0_and_d12_or_c2 == 1 && dst_pow2 <= 81.0f) {
	hpc1 += myprt->prt_point[p].hpp *
	  (1.0f - (3.5f / 81.0f * dst_pow2 -
		   4.5f / 81.0f / 81.0f * dst_pow4 +
		   2.5f / 81.0f / 81.0f / 81.0f * dst_pow4 * dst_pow2 -
		   0.5f / 81.0f / 81.0f / 81.0f / 81.0f * dst_pow4 * dst_pow4));
      }
      

      /* L-J potential */
      const float p1 = enepara_dc->pa[lig_t][prt_t][0] / (dst_pow4 * dst_pow4 * dst);
      const float p2 = enepara_dc->pa[lig_t][prt_t][1] / (dst_pow4 * dst_pow2);
      const float p4 = p1 * enepara_lj0_dc * (1.0f + enepara_lj1_dc * dst_pow2) + 1.0f;
      evdw += (p1 - p2) / p4;




      /* electrostatic potential */
      const float s1 = enepara_el1_dc * dst;
      float g1;
      if (s1 < 1)
	g1 = enepara_el0_dc + enepara_a1_dc * s1 * s1 + enepara_b1_dc * s1 * s1 * s1;
      else
	g1 = 1.0f / s1;
      eele += mylig->lig_point[l].c * myprt->prt_point[p].ele * g1;

      
      /* contact potential */
      const float dst_minus_pmf0 = dst - enepara_dc->pmf[lig_t][prt_t][0];

      epmf +=
	enepara_dc->pmf[lig_t][prt_t][1] /
	(1.0f + expf ((-0.5f * dst + 6.0f) * dst_minus_pmf0));



      /* pocket-specific potential */
      // the senmatics do not match with the original program:
      // if (found psp[][])
      //   accumulate to epsp;
      // else
      //   do nothing
      if (myprt->prt_point[p].c == 2 && dst_minus_pmf0 <= 0) {
	const int i1 = myprt->prt_point[p].seq3r;
	epsp += psp_dc->psp[lig_t][i1]; // sparse matrix
      }


      /* hydrogen bond potential */
      const float hdb0 = enepara_dc->hdb[lig_t][prt_t][0];
      if (hdb0 > 0.1f) {
	const float hdb1 = enepara_dc->hdb[lig_t][prt_t][1];
	const float hdb3 = (dst - hdb0) * hdb1;
	ehdb += hdb1 * expf (-0.5f * hdb3 * hdb3);
      }


    } // prt loop



    /* hydrophobic restraits*/
    // transpose may help improve the performance
    const float hpc2 = (hpc1 - enepara_dc->hpl[lig_t][0]) / enepara_dc->hpl[lig_t][1];
    ehpc += 0.5f * hpc2 * hpc2 - enepara_dc->hpl[lig_t][2];

  } // lig loop







#if 1

  /* kde potential */


/* icc -mmic -vec-report3: loop was not vectorized: not inner loop */
  //#pragma omp parallel for reduction(+:ekde)
  // lig loop, ~30
  for (int l = 0; l < lna_dc; ++l) {
    float kde_val = 0.0f;
    float kde_sz = 0.0f;
    const int lig_t = mylig->lig_point[l].t;

/* icc -mmic -vec-report3: LOOP WAS VECTORIZED */
    // kde loop, ~400
    for (int k = 0; k < pnk_dc; ++k) {
      if (lig_t == kde_dc->kde_point[k].t) {
	const float dx = mylig->coord_new.x[l] - kde_dc->kde_point[k].x;
	const float dy = mylig->coord_new.y[l] - kde_dc->kde_point[k].y;
	const float dz = mylig->coord_new.z[l] - kde_dc->kde_point[k].z;
	const float kde_dst_pow2 = dx * dx + dy * dy + dz * dz;
	kde_val += expf (enepara_kde2_dc * kde_dst_pow2);
	kde_sz += 1.0f;
      }
    } // kde loop

    if (kde_sz != 0.0f)
      ekde += (kde_val / kde_sz);

  } // lig loop

  ekde = ekde / enepara_kde3_dc;


#endif








#if 1

  /* position restraints */

/* icc -mmic -vec-report3: loop was not vectorized: not inner loop */
  //#pragma omp parallel for reduction(+:elhm)
  // lhm loop, ~11
  for (int m = 0; m < n_pos_dc; ++m) {
    float lhm_val = 0.0f;
    float lhm_sz = 0.0f;

/* icc -mmic -vec-report3: LOOP WAS VECTORIZED */
/* icc -mmic -vec-report3: PEEL LOOP WAS VECTORIZED */
    // lig loop, ~30
    for (int l = 0; l < lna_dc; ++l) {
      const int lig_n = mylig->lig_point[l].n + 1;
      if (mcs_dc[m].mcs_point[lig_n].x != MCS_INVALID_COORD) {
	const float dx = mylig->coord_new.x[l] - mcs_dc[m].mcs_point[lig_n].x;
	const float dy = mylig->coord_new.y[l] - mcs_dc[m].mcs_point[lig_n].y;
	const float dz = mylig->coord_new.z[l] - mcs_dc[m].mcs_point[lig_n].z;
	lhm_val += dx * dx + dy * dy + dz * dz;
	lhm_sz += 1.0f;
      }
    } // lig loop

    elhm += mcs_dc[m].tcc * sqrtf (lhm_val / lhm_sz);

  } // lhm loop


#endif


  // energy edst e[8]
  const float dx = mylig->coord_new.center[0] - myprt->pocket_center[0];
  const float dy = mylig->coord_new.center[1] - myprt->pocket_center[1];
  const float dz = mylig->coord_new.center[2] - myprt->pocket_center[2];
  float edst = sqrtf (dx * dx + dy * dy + dz * dz);



  Energy e;
  e.e[0] = evdw / lna_dc;
  e.e[1] = eele / lna_dc;
  e.e[2] = epmf / lna_dc;
  e.e[3] = ehpc / lna_dc;
  e.e[4] = ehdb / lna_dc / sqrtf (2.0f * PI) * -1.0f;
  e.e[5] = edst;
  e.e[6] = epsp / lna_dc;
  e.e[7] = ekde / lna_dc;
  e.e[8] = logf (elhm / n_pos_dc);
   
  // calculate normalized energy
  for (int i = 0; i < MAXWEI - 1; ++i)
    e.e[i] = enepara_dc->a_para[i] * e.e[i] + enepara_dc->b_para[i];

  // calculate the total energy using linear combination
  float etotal = 0.0f;
  for (int i = 0; i < MAXWEI - 1; ++i)
    etotal +=  enepara_dc->w[i] * e.e[i];
  e.e[MAXWEI - 1] = etotal;

  mylig->energy_new = e;

}


