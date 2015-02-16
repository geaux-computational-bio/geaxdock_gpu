#include <cmath>
#include <cstdio>

#include "dock.h"


extern int n_pos;


void
CalcEnergy (Ligand *lig, Protein *prt, Psp *psp, Kde *kde, Mcs *mcs, EnePara *enepara)
{
  Ligand *mylig = &lig[0];
  Protein *myprt = &prt[0];
  Lig_Coord *myligcoord = &mylig->coord[mylig->track];


  // Energies
  float evdw = 0.0f;
  float eele = 0.0f;
  float epmf = 0.0f;
  float ehpc = 0.0f;
  float ehdb = 0.0f;
  float epsp = 0.0f;
  float ekde = 0.0f;
  float elhm = 0.0f;



  for (int i = 0; i < mylig->lna; ++i) {
    int lig_t = mylig->t[i];
    float hpc1 = 0.0f;


    for (int j = 0; j < myprt->pnp; ++j) {
      int prt_t = myprt->t[j];

      float dst_pow2 =
	powf (myligcoord->x[i] - myprt->x[j], 2.0f) +
	powf (myligcoord->y[i] - myprt->y[j], 2.0f) +
	powf (myligcoord->z[i] - myprt->z[j], 2.0f);

      float dst = sqrtf (dst_pow2);



      /* L-J potential */
      float p1 = enepara->p1a[lig_t][prt_t] / powf (dst, 9.0f);
      float p2 = enepara->p2a[lig_t][prt_t] / powf (dst_pow2, 3.0f);
      float p4 = p1 * enepara->lj0 * (1.0f + enepara->lj1 * dst_pow2) + 1.0f;
      evdw += (p1 - p2) / p4;




      /* electrostatic potential */

      float s1 = enepara->el1 * dst;
      float g1;

      if (s1 < 1)
	g1 = enepara->el0 +
	  enepara->a1 * powf (s1, 2.0f) +
	  enepara->b1 * powf (s1, 3.0f);
      else
	g1 = 1.0f / s1;

      int e3_idx = prt_t == 0 ? myprt->d[j] + 30 : prt_t;
      float e3 = enepara->ele[e3_idx];

      // energy eele e[1]
      eele += mylig->c[i] * e3 * g1;




      /* contact potential */
      // energy epmf e[2]
      epmf +=
	enepara->pmf1[lig_t][prt_t] *
	1.0f / (1.0f + expf ((-0.5f * dst + 6.0f) *
			   (dst - enepara->pmf0[lig_t][prt_t])));



      /* hydrogen bond potential */
      float tmp0 = enepara->hdb0[lig_t][prt_t];
      float tmp1 = enepara->hdb1[lig_t][prt_t];
      if (tmp0 > 0.1f)
	// energy ehdb e[4]
	ehdb += (-1.0f / (tmp1 * sqrtf (2.0f * PI))) *
	  expf (-0.5f * powf ((dst - tmp0) / tmp1, 2.0f));





      // not exactly match the original program
      /* pocket-specific potential */
      if (myprt->c[j] == 2 && dst <= enepara->pmf0[lig_t][prt_t]) {
	int pspidx2 = myprt->seq3r[j];
	// energy epsp e[5]
	epsp += psp->psp[lig_t][pspidx2]; // sparse matrix

	//printf ("access psp \t%2d\t%2d\n", lig_t, pspidx2);
      }

      //printf ("unconditional access psp \t%2d\t%2d", lig_t, myprt->seq3r[j]);
      //printf ("\t prt_c = %2d\n", myprt->c[j]);






      /* hydrophobic potential */
      if (myprt->c0_and_d12_or_c2[j] && dst <= 9.0f) {

#if 0
	hpc1 += enepara->hpp[myprt->d[j]] *
	  (1.0f - 0.5f * (7.0f * powf (dst / 9.0f, 2.0f) -
			  9.0f * powf (dst / 9.0f, 4.0f) +
			  5.0f * powf (dst / 9.0f, 6.0f) -
			  powf (dst / 9.0f, 8.0f)));
#endif

	float dst9_pow2 = powf (dst / 9.0f, 2.0f);
	hpc1 += enepara->hpp[myprt->d[j]] *
	  (1.0f - 0.5f * (7.0f * dst9_pow2 -
			  9.0f * powf (dst9_pow2, 2.0f) +
			  5.0f * powf (dst9_pow2, 3.0f) -
			  powf (dst9_pow2, 4.0f)));


      }






    }				// protein loop, j loop



    /* hydrophobic restraints */
    // energy ehpc e[3]
    ehpc +=
      0.5f * powf ((hpc1 - enepara->hpl0[lig_t]) / enepara->hpl1[lig_t], 2.0f) -
      log (1.0f / (enepara->hpl1[lig_t] * sqrtf (2.0f * PI)));

  }				// ligand loop, i loop








  /* kde potential */
  for (int i = 0; i < mylig->lna; ++i) {
    float lig_x = myligcoord->x[i];
    float lig_y = myligcoord->y[i];
    float lig_z = myligcoord->z[i];
    int lig_t = mylig->t[i];

    float kde_sum = 0.0f;
    int kde_size = 0;

    for (int j = 0; j < kde->pnk; ++j) {
      if (lig_t == kde->t[j]) {
	float tmp1 = enepara->kde2;
	float tmp2 = tmp1 * sqrtf (2.0f * PI);
	float tmp3 =
	  powf ((lig_x - kde->x[j]) / tmp1, 2.0f) +
	  powf ((lig_y - kde->y[j]) / tmp1, 2.0f) +
	  powf ((lig_z - kde->z[j]) / tmp1, 2.0f);
	kde_sum += expf (-0.5f * tmp3) / (tmp2 * tmp2 * tmp2);
	kde_size++;

	/*
	kde_sum +=
	  ((expf (-0.5f * powf ((lig_x - kde->x[j]) / tmp1, 2.0f)) / tmp2) *
	   (expf (-0.5f * powf ((lig_y - kde->y[j]) / tmp1, 2.0f)) / tmp2) *
	   (expf (-0.5f * powf ((lig_z - kde->z[j]) / tmp1, 2.0f)) / tmp2));
	*/

      }
    }
				// kde loop, j loop
    if (kde_size != 0)
      // energy ekde e[6]
      ekde += (kde_sum / kde_size);

  }				// ligand loop, i loop






#if 0
  /* position restraints */

  for (int mcs_idx = 0; mcs_idx < n_pos; ++mcs_idx) {

    float rmsd1_sum = 0.0f;
    int rmsd1_size = 0;

    // ~30
    for (int lig_idx = 0; lig_idx < mylig->lna; ++lig_idx) {
      int lig_n = mylig->n[lig_idx] + 1;
      float mcs_x = mcs->x[mcs_idx][lig_n];

      if (mcs_x != MCS_INVALID_COORD) {
	float lig_x = myligcoord->x[lig_idx];
	float lig_y = myligcoord->y[lig_idx];
	float lig_z = myligcoord->z[lig_idx];
        float mcs_y = mcs->y[mcs_idx][lig_n];
        float mcs_z = mcs->z[mcs_idx][lig_n];
	rmsd1_sum +=
	  (powf (lig_x - mcs_x, 2.0f) +
	   powf (lig_y - mcs_y, 2.0f) +
	   powf (lig_z - mcs_z, 2.0f));
	rmsd1_size += 1;
      }

    }

    // energy elhm e[7]
    elhm += mcs->tcc[mcs_idx] * sqrtf (rmsd1_sum / rmsd1_size);
  }



  if (n_pos != 0)
    elhm = logf (elhm / n_pos);
#endif






  // energy edst e[8]
  float edst =
    sqrtf (powf (myligcoord->center[0] - mylig->pocket_center[0], 2.0f) +
	  powf (myligcoord->center[1] - mylig->pocket_center[1], 2.0f) +
	  powf (myligcoord->center[2] - mylig->pocket_center[2], 2.0f));





  Energy energy;
  energy.e[0] = evdw;
  energy.e[1] = eele;
  energy.e[2] = epmf;
  energy.e[3] = ehpc;
  energy.e[4] = ehdb;
  energy.e[5] = epsp;
  energy.e[6] = ekde;
  energy.e[7] = elhm;
  energy.e[8] = edst;




  float lna = (float) mylig->lna;
  for (int i = 0; i <= 6; ++i)	// not for elhm, edst
    energy.e[i] /= lna;


  float etotal = 0.0f;
  for (int i = 0; i < MAXWEI; ++i)
    etotal += enepara->w[i] * energy.e[i];
  energy.total = etotal;


  // save energy
  mylig->energy[mylig->track] = energy;
}


