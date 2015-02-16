#include <cstdio>

#include "dock.h"




void
MonteCarlo_Init (Ligand *lig, Protein *prt, Psp *psp, Kde *kde, Mcs *mcs, EnePara *enepara)
{
  int lig_idx = 0;
  Ligand *mylig = &lig[lig_idx];

  mylig->track = mylig->track ^ 1;
  CalcEnergy (lig, prt, psp, kde, mcs, enepara);

  PrintEnergy2 (&mylig->energy[mylig->track], 999, 3);

  mylig->track = mylig->track ^ 1;
}





void
MonteCarlo (Ligand *lig, Protein *prt, Psp *psp, Kde *kde, Mcs *mcs, EnePara *enepara, McPara *mcpara)
{
  //const float minus_beta = mcpara->minus_beta;
  const int step = mcpara->step;

  int lig_idx = 0;
  Ligand *mylig = &lig[lig_idx];

  for (int i = 0; i < step; ++i) {
    Move (lig, mcpara->t, mcpara->r);
    CalcEnergy (lig, prt, psp, kde, mcs, enepara);
    PrintEnergy2 (&mylig->energy[mylig->track], i, 2);
    Accept (lig, mcpara);

/*
    printf ("%08d \t track = %d \t new_energy = %f \t old_energy = %f delta_energy = %f \t is_accept = %d\n",
	    i,
	    track,
	    mylig->energy[track].total,
	    mylig->energy[!track].total,
	    delta_energy,
	    is_accept);
*/
  }
}


