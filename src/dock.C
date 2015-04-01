#include <iostream>
#include <cstdlib>
#include <fstream>
#include <ctime>
#include <string>
#include <cstdio>

#include "size.h"
#include "dock.h"
#include "run.h"
#include "util.h"
#include "toggle.h"
#include "load.h"
#include "stats.h"
#include "post_mc.h"


using namespace std;

int
main (int argc, char **argv)
{
  Banner ();

  srand (time (0));

  McPara *mcpara = new McPara;
  McLog *mclog = new McLog;
  ExchgPara *exchgpara = new ExchgPara;
  InputFiles *inputfiles = new InputFiles[1];
  ParseArguments (argc, argv, mcpara, exchgpara, inputfiles);

  // load into preliminary data structures
  Ligand0 *lig0 = new Ligand0[MAXEN2];
  Protein0 *prt0 = new Protein0[MAXEN1];
  Psp0 *psp0 = new Psp0;
  Kde0 *kde0 = new Kde0;
  Mcs0 *mcs0 = new Mcs0[MAXPOS];
  EnePara0 *enepara0 = new EnePara0;

  loadLigand (inputfiles, lig0);
  loadProtein (&inputfiles->prt_file, prt0);
  loadLHM (&inputfiles->lhm_file, psp0, kde0, mcs0);
  loadEnePara (&inputfiles->enepara_file, enepara0);
  

  // sizes
  ComplexSize complexsize;
  complexsize.n_prt = inputfiles->prt_file.conf_total;	// number of protein conf
  complexsize.n_tmp = exchgpara->num_temp;	// number of temperature
  complexsize.n_lig = inputfiles->lig_file.conf_total;	// number of ligand conf
  complexsize.n_rep = complexsize.n_lig * complexsize.n_prt * complexsize.n_tmp;
  complexsize.lna = inputfiles->lig_file.lna;
  complexsize.pnp = inputfiles->prt_file.pnp;
  complexsize.pnk = kde0->pnk;
  complexsize.pos = inputfiles->lhm_file.pos;	// number of MCS positions


  // data structure optimizations 
  Ligand *lig = new Ligand[complexsize.n_rep];
  Protein *prt = new Protein[complexsize.n_prt];
  Psp *psp = new Psp;
  Kde *kde = new Kde;
  Mcs *mcs = new Mcs[complexsize.pos];
  EnePara *enepara = new EnePara;
  Temp *temp = new Temp[complexsize.n_tmp];
  Replica *replica = new Replica[complexsize.n_rep];

  OptimizeLigand (lig0, lig, complexsize);
  OptimizeProtein (prt0, prt, enepara0, lig0, complexsize);
  OptimizePsp (psp0, psp, lig, prt);
  OptimizeKde (kde0, kde);
  OptimizeMcs (mcs0, mcs, complexsize);
  OptimizeEnepara (enepara0, enepara);

  delete[]lig0;
  delete[]prt0;
  delete[]psp0;
  delete[]kde0;
  delete[]mcs0;
  delete[]enepara0;


  // initialize system
  InitLigCoord (lig, complexsize);
  SetTemperature (temp, exchgpara);
  // SetTemperature (temp, mcpara, complexsize);
  SetReplica (replica, lig, complexsize);
  SetMcLog (mclog);

  // debug
  //PrintDataSize (lig, prt, psp, kde, mcs, enepara);
  //PrintLigand (lig);
  //PrintProtein (prt);


  map < int, vector < LigRecordSingleStep > > multi_reps_records;
  printf ("Start docking\n");
  Run (lig, prt, psp, kde, mcs, enepara, temp, replica, mcpara, mclog,
       multi_reps_records, complexsize);

  // post_mc(multi_reps_records, lig, prt, enepara, mcpara);
  opt_ff(multi_reps_records, lig, complexsize.n_lig, prt, enepara, mcpara);

  PrintSummary (inputfiles, mcpara, temp, mclog, &complexsize);

  // clean up

  delete[]mcpara;
  delete[]mclog;
  delete[]inputfiles;
  delete[]lig;
  delete[]prt;
  delete[]psp;
  delete[]kde;
  delete[]mcs;
  delete[]enepara;
  delete[]temp;
  delete[]replica;
  delete[]exchgpara;

  return 0;
}
