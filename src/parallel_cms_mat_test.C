#include <cstdio>
#include <iostream>
#include <vector>
#include <random>

#include <omp.h>

#include "load.h"
#include "dock.h"
#include "size.h"
#include "util.h"

#include "gtest/gtest.h"
#include "gtest/internal/gtest-internal.h"

using namespace std;

void
printMatrix(double** mat, int n)
{
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j)
      {
        printf("%.3f ", mat[i][j]);
      }
    cout << endl;
  }
}

TEST(OpenMp, Simple)
{
#pragma omp parallel num_threads(3)
  {
    int i = omp_get_thread_num();
    cout << "open mp" << i << endl;
  }
}

TEST(load, complex)
{
  McPara *mcpara = new McPara;
  McLog *mclog = new McLog;
  ExchgPara *exchgpara = new ExchgPara;
  InputFiles *inputfiles = new InputFiles[1];
  
  inputfiles->lig_file.path = "../data/1robA1/1robA1.sdf";
  inputfiles->prt_file.path = "../data/1robA1/1robA.pdb";
  inputfiles->lhm_file.path = "../data/1robA1/1robA1-0.8.ff";
  inputfiles->enepara_file.path = "../data/parameters/paras";
  
  exchgpara->num_temp = 1;

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
  SetReplica (replica, lig, complexsize);
  
  int tot_steps = 300;
  vector < LigRecordSingleStep > steps(tot_steps);

  float HI = 0.1, LO = -0.1;
  for (auto iter = steps.begin(); iter != steps.end(); ++iter) {
    float r3 = LO + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(HI-LO)));
    for (int i = 0; i < 6; ++i)
      (*iter).movematrix[i] = r3;
  }
  
  
  int tot = steps.size();
  


  double** dis_mat = AllocSquareMatrix(tot);
  double** p_dis_mat = AllocSquareMatrix(tot);

  double starting_time = get_wall_time();
  // Serial version
  GenCmsSimiMat(steps, lig, prt, enepara, dis_mat);
  double first_now = get_wall_time();

  cout << "Serial version takes: " << first_now - starting_time << endl;
  
  // OpenMp version
  int total_threads = 4;
  int n_lig = complexsize.n_lig;
  ParallelGenCmsSimiMat(steps, lig, n_lig, prt, enepara, p_dis_mat);
  double second_now = get_wall_time();
  cout << "OpenMp version takes: " << second_now - first_now << endl;

  
  // compare results
  for (int i = 0; i < tot; ++i)
    for (int j = 0; j < tot; ++j)
      {
        assert((dis_mat[i][j] - p_dis_mat[i][j]) < 0.001);
      }
  
  
  FreeSquareMatrix(dis_mat, tot);
  FreeSquareMatrix(p_dis_mat, tot);

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
}
