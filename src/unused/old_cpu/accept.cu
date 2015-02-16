#include <cmath>
#include <cstdlib>
#include <cstdio>

#include "dock.h"

#include <cuda.h>
#include "gpu.cuh"



void
Accept (Ligand *lig, McPara* mcpara)
{
  int lig_idx = 0;
  Ligand *mylig = &lig[lig_idx];


  //int track = mylig->track;
  //float delta_energy = (mylig->energy[track].total) - (mylig->energy[!track].total);
  int is_accept = 1;
//    int is_accept =
//      (float) rand () / RAND_MAX < expf (delta_energy * minus_beta);
  mylig->track ^= is_accept;
  mcpara->ar += is_accept;
}




