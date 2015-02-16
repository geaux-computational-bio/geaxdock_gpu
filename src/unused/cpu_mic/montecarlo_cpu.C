/*
#include <cstdio>

#include "dock.h"

#include <cuda.h>
#include "gpu.cuh"
*/



void
MonteCarlo_Init_d (const int rep_begin, const int rep_end)
{
  for (int myreplica = rep_begin; myreplica <= rep_end; ++myreplica) {
    Ligand *mylig = &lig_dc[replica_dc[myreplica].idx_rep];
    const Protein *myprt = &prt_dc[replica_dc[myreplica].idx_prt];
    CalcEnergy_d (mylig, myprt);
    mylig->energy_old = mylig->energy_new;
  }

  InitAcs_d (); // set acceptance counters to zeros

}



void
MonteCarlo_d (const int rep_begin, const int rep_end, const int s1, const int s2)
{

//#pragma offload target(mic)
{
#pragma omp parallel for
  for (int myreplica = rep_begin; myreplica <= rep_end; ++myreplica) {
    Ligand *mylig = &lig_dc[replica_dc[myreplica].idx_rep];
    const Protein *myprt = &prt_dc[replica_dc[myreplica].idx_prt];
    const float mybeta = temp_dc[replica_dc[myreplica].idx_tmp].minus_beta;
    // printf("mybeta: %f\n", mybeta);
    
    for (int s3 = 0; s3 < steps_per_exchange_dc; ++s3) {

#if IS_OUTPUT == 1
      // record old status
      RecordLigand_d (s1, s2 + s3, myreplica, rep_begin, mylig);
#endif
      Move_d (mylig);
      CalcEnergy_d (mylig, myprt);
      Accept_d (mylig, mybeta, myreplica);
    }

    etotal_dc[myreplica] = mylig->energy_old.e[MAXWEI - 1];
    for (int i = 0; i < 6; ++i)
      ligmovevector_dc[myreplica].ele[i] = mylig->movematrix_old[i];
  }
}

}
