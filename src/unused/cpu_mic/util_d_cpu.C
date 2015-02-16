/*
#include <cstdlib>
#include <cstdio>

#include "dock.h"
#include "gpu.cuh"

#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
*/





void
InitAcs_d ()
{
  for (int i = 0; i < MAXREP; ++i) {
    acs_dc[i] = 0;
  }
}

void
InitLigRecord_d (const int myreplica, const int rep_begin)
{
  for (int s2s3 = 0; s2s3 < steps_per_dump_dc; ++s2s3) {
    LigRecordSingleStep *myrecord = &ligrecord_dc[myreplica - rep_begin].step[s2s3];
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



void
RecordLigand_d (const int s1, const int s2s3, const int myreplica, const int rep_begin, const Ligand * mylig)
{
  LigRecordSingleStep *myrecord = &ligrecord_dc[myreplica - rep_begin].step[s2s3];
  myrecord->replica = replica_dc[myreplica];
  myrecord->energy = mylig->energy_old;
  for (int i = 0; i < 6; ++i)
    myrecord->movematrix[i] = mylig->movematrix_old[i];
  myrecord->step = s1 + s2s3;
}



float
MyRand_d ()
{

  float randdd; 

  if (is_random_dc == 0) {
     randdd = 20.0f;
    //randdd = 0.0f;
  }
  else {
    randdd = (float) rand() / (float) RAND_MAX;
  }

  // printf("%f\n", randdd);

  return randdd;
}


