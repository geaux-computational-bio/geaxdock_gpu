
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
CombineEnergy_d (const int bidx, Energy * e)
{
  if (bidx == 0) {

  // calculate the total energy using linear combination
#if IS_LINEAR == 1
    float etotal = 0.0f;
    for (int i = 0; i < MAXWEI - 1; ++i) {
      float * ener = &e->e[i];
      *ener = (*ener) * enepara_dc->w[i];
      etotal += *ener;
    }
    e->e[MAXWEI - 1] = etotal;
#endif


  // calculate the total energy using Bayes' formula
#if IS_BAYE == 1
#include "distribution.h"
    float eh[MAXWEI], el[MAXWEI]; //conditional prob belonging to high decoy

    // 0 - vdw
    eh[0] = NormPdf(e->e[0], VDW_NORM_HIGH_LOC, VDW_NORM_HIGH_SCALE);
    el[0] = NormPdf(e->e[0], VDW_NORM_LOW_LOC, VDW_NORM_LOW_SCALE);
    // 1 - ele
    eh[1] = CauchyPdf(e->e[1], ELE_CAUCHY_HIGH_LOC, ELE_CAUCHY_HIGH_SCALE);
    el[1] = CauchyPdf(e->e[1], ELE_CAUCHY_LOW_LOC, ELE_CAUCHY_LOW_SCALE);
    // 2 - pmf
    eh[2] = LogisticPdf(e->e[2], PMF_LOGISTIC_HIGH_LOC, PMF_LOGISTIC_HIGH_SCALE);
    el[2] = LogisticPdf(e->e[2], PMF_LOGISTIC_LOW_LOC, PMF_LOGISTIC_LOW_SCALE);
    // 3 - psp
    eh[3] = LogisticPdf(e->e[3], PSP_LOGISTIC_HIGH_LOC, PSP_LOGISTIC_HIGH_SCALE);
    el[3] = LogisticPdf(e->e[3], PSP_LAPLACE_LOW_LOC, PSP_LAPLACE_LOW_SCALE);
    // 4 - hdb
    eh[4] = NormPdf(e->e[4], HDB_NORM_HIGH_LOC, HDB_NORM_HIGH_SCALE);
    el[4] = NormPdf(e->e[4], HDB_LOGISTIC_LOW_LOC, HDB_LOGISTIC_LOW_SCALE);
    // 5 - hpc
    eh[5] = WaldPdf(e->e[5], HPC_WALD_HIGH_LOC, HPC_WALD_HIGH_SCALE);
    el[5] = WaldPdf(e->e[5], HPC_WALD_LOW_LOC, HPC_WALD_LOW_SCALE);
    // 6 - kde
    eh[6] = WaldPdf(e->e[6], KDE_WALD_HIGH_LOC, KDE_WALD_HIGH_SCALE);
    el[6] = WaldPdf(e->e[6], KDE_WALD_LOW_LOC, KDE_WALD_LOW_SCALE);
    // 7 - lhm
    eh[7] = LogisticPdf(e->e[7], LHM_LOGISTIC_HIGH_LOC, LHM_LOGISTIC_HIGH_SCALE);
    el[7] = LogisticPdf(e->e[7], LHM_LOGISTIC_LOW_LOC, LHM_LOGISTIC_LOW_SCALE);
    // 8 - dst
    eh[8] = LogisticPdf(e->e[8], DST_LOGISTIC_HIGH_LOC, DST_LOGISTIC_HIGH_SCALE);
    el[8] = LogisticPdf(e->e[8], DST_LOGISTIC_LOW_LOC, DST_LOGISTIC_LOW_SCALE);
    
    // calculate conditional prob
    float prob_h = 0.0f, prob_l = 0.0f;
    for (int i = 0; i <= 8; ++i) {
      prob_h += log10f(eh[i]);
      prob_l += log10f(el[i]);
    }
    e->e[MAXWEI - 1] = prob_l - prob_h;
#endif

  }

}


