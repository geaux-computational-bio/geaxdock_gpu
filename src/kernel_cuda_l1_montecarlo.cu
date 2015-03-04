/*
#include <cstdio>

#include "dock.h"

#include <cuda.h>
#include "gpu.cuh"
*/



__global__ void
MonteCarlo_Init_d (const int rep_begin, const int rep_end)
{
  const int bidx = blockDim.x * threadIdx.y + threadIdx.x;      // within a TB

  for (int offset = rep_begin; offset <= rep_end; offset += GD) {
    const int myreplica = offset + blockIdx.x;

    if (myreplica <= rep_end) {
      Ligand *mylig = &lig_dc[replica_dc[myreplica].idx_rep];
      const Protein *myprt = &prt_dc[replica_dc[myreplica].idx_prt];

      if (myreplica == 0)
        InitRefMatrix_d (bidx, mylig, myprt);
        // mcc ref matrix generated from the first replica

#if IS_AWAY == 1
      Move_d (bidx, mylig, 44.5f);
#endif 
      CalcEnergy_d (bidx, mylig, myprt);
	   
#if IS_CALCU_RMSD == 1
      CalcRmsd_d (bidx, mylig);
#endif

#if IS_CALCU_MCC == 1
      CalcMcc_d (bidx, mylig, myprt);
#endif

      if (bidx < MAXWEI)
	mylig->energy_old.e[bidx] = mylig->energy_new.e[bidx];

      if (bidx == 0)
	mylig->energy_old.cms = mylig->energy_new.cms;
      
#if IS_AWAY
      // force to accept, set mybeta to be zero
      Accept_d (bidx, mylig, 0.000000f, myreplica);
#endif

      if (bidx == 0)
	mylig->is_move_accepted = 1;

#if IS_OUTPUT == 1
	// record old status
      RecordLigand_d (bidx, 0, 0, myreplica, rep_begin, mylig);
#endif


    }
  }

  InitAcs_d (bidx); // set acceptance counters to zeros

}



__global__ void
MonteCarlo_d (const int rep_begin, const int rep_end, const int s1, const int s2)
{
  const int bidx = blockDim.x * threadIdx.y + threadIdx.x;      // within a TB

  //__shared__ LigCoord myligcoord[1];

  for (int offset = rep_begin; offset <= rep_end; offset += GD) {
    const int myreplica = offset + blockIdx.x;

    /*
    if (bidx == 0) {
      printf ("%3d : %3d\n", myreplica, replica_dc[myreplica].idx_rep);
    }
    */

    if (myreplica <= rep_end) {
      Ligand *mylig = &lig_dc[replica_dc[myreplica].idx_rep];
      const Protein *myprt = &prt_dc[replica_dc[myreplica].idx_prt];
      const float mybeta = temp_dc[replica_dc[myreplica].idx_tmp].minus_beta;
      // printf("mybeta: %f\n", mybeta);


      for (int s3 = 0; s3 < steps_per_exchange_dc; ++s3) {
	
#if IS_CONTROL_MOVE == 1
	Move_d (bidx, mylig, 2.0f);
#else
	Move_d (bidx, mylig, 2.0f * MyRand_d() - 1.0f);
#endif

#if IS_CALCU_RMSD == 1
        CalcRmsd_d (bidx, mylig);
#endif

#if IS_CALCU_MCC == 1
	CalcMcc_d (bidx, mylig, myprt);
#endif 

	CalcEnergy_d (bidx, mylig, myprt);
	Accept_d (bidx, mylig, mybeta, myreplica);

#if IS_OUTPUT == 1
	// record old status
	RecordLigand_d (bidx, s1, s2 + s3, myreplica, rep_begin, mylig);
#endif
      }

      if (bidx == 0) {
	etotal_dc[myreplica] = mylig->energy_old.e[MAXWEI - 1];
	for (int i = 0; i < 6; ++i)
	  ligmovevector_dc[myreplica].ele[i] = mylig->movematrix_old[i];
      }
    }
  }

}

