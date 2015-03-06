//#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <assert.h>

#include "dock.h"
#include "size.h"
#include "toggle.h"
#include "hdf5io.h"

#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>

#include "dock.h"
#include "toggle.h"
#include "util.h"
#include "kernel_cuda.cuh"


#include <yeah/cudautil.h>
#include <yeah/timing.h>



// inline everything in order to use constant memory
#include "kernel_cuda.cu"



void
Run (const Ligand * lig,
     const Protein * prt,
     const Psp * psp,
     const Kde * kde,
     const Mcs * mcs,
     const EnePara * enepara,
     const Temp * temp,
     const Replica * replica,
     const McPara * mcpara,
     McLog * mclog,
     const ComplexSize complexsize)
{
  //Parameter para;


  // sizes
  const int n_lig = complexsize.n_lig;
  const int n_prt = complexsize.n_prt;
  const int n_tmp = complexsize.n_tmp;
  const int n_rep = complexsize.n_rep;



  // set up kernel parameters
  const int dim_grid = GD;
  const dim3 dim_block (BDx, BDy, 1);



  // initilize random sequence
  srand (time (NULL));
  curandState *curandstate_d[NGPU];
  size_t curandstate_sz = sizeof (curandState) * TperB * GD;

  for (int i = 0; i < NGPU; ++i) {
    cudaSetDevice (i);
    int myrand = rand () + i;
    CUDAMEMCPYTOSYMBOL (seed_dc, &myrand, int);

    CUDAMALLOC (curandstate_d[i], curandstate_sz, curandState *);
    CUDAMEMCPYTOSYMBOL (curandstate_dc, &curandstate_d[i], curandState *);
    CUDAKERNELSYNC (InitCurand_d, dim_grid, dim_block);
  }




  //int ngpu = 0;
  //cudaGetDeviceCount (&ngpu);

  for (int i = 0; i < NGPU; ++i) {
    cudaSetDevice (i);
    cudaFuncSetCacheConfig (MonteCarlo_Init_d, cudaFuncCachePreferShared);
    cudaFuncSetCacheConfig (MonteCarlo_d, cudaFuncCachePreferShared);
    cudaFuncSetCacheConfig (ExchangeReplicas_d, cudaFuncCachePreferShared);
  }


  // calculate the upper/lower bound for multi-GPU data decomposition
  const int n_rep_per_gpu_max = (int) ceilf ((float) (n_prt * n_tmp) / NGPU) * n_lig;
  const int n_active_gpu = (int) ceilf ((float) n_rep / n_rep_per_gpu_max);
  if (n_active_gpu < NGPU) {
    fprintf (stderr, "error: n_active_gpu < NGPU");
    exit (9999);
  }

  int rep_begin[NGPU], rep_end[NGPU], n_rep_per_gpu[NGPU];
  for (int i = 0; i < NGPU; ++i) {
    rep_begin[i] = n_rep_per_gpu_max * i;
    rep_end[i] = minimal_int (rep_begin[i] + n_rep_per_gpu_max - 1, n_rep);
    n_rep_per_gpu[i] = rep_end[i] - rep_begin[i] + 1;
  }




  // GPU read only scalars
  for (int i = 0; i < NGPU; ++i) {
    cudaSetDevice (i);
    CUDAMEMCPYTOSYMBOL (steps_total_dc, &mcpara->steps_total, int);
    CUDAMEMCPYTOSYMBOL (steps_per_dump_dc, &mcpara->steps_per_dump, int);
    CUDAMEMCPYTOSYMBOL (steps_per_exchange_dc, &mcpara->steps_per_exchange, int);

    CUDAMEMCPYTOSYMBOL (enepara_lj0_dc, &enepara->lj0, float);
    CUDAMEMCPYTOSYMBOL (enepara_lj1_dc, &enepara->lj1, float);
    CUDAMEMCPYTOSYMBOL (enepara_el0_dc, &enepara->el0, float);
    CUDAMEMCPYTOSYMBOL (enepara_el1_dc, &enepara->el1, float);
    CUDAMEMCPYTOSYMBOL (enepara_a1_dc, &enepara->a1, float);
    CUDAMEMCPYTOSYMBOL (enepara_b1_dc, &enepara->b1, float);
    CUDAMEMCPYTOSYMBOL (enepara_kde2_dc, &enepara->kde2, float);
    CUDAMEMCPYTOSYMBOL (enepara_kde3_dc, &enepara->kde3, float);

    CUDAMEMCPYTOSYMBOL (n_lig_dc, &n_lig, int);
    CUDAMEMCPYTOSYMBOL (n_prt_dc, &n_prt, int);
    CUDAMEMCPYTOSYMBOL (n_tmp_dc, &n_tmp, int);
    CUDAMEMCPYTOSYMBOL (n_rep_dc, &n_rep, int);
    CUDAMEMCPYTOSYMBOL (lna_dc, &complexsize.lna, int);
    CUDAMEMCPYTOSYMBOL (pnp_dc, &complexsize.pnp, int);
    CUDAMEMCPYTOSYMBOL (pnk_dc, &complexsize.pnk, int);
    CUDAMEMCPYTOSYMBOL (pos_dc, &complexsize.pos, int);
  }



  // GPU read only arrays
  const size_t prt_sz = sizeof (Protein) * n_prt;
  const size_t psp_sz = sizeof (Psp);
  const size_t kde_sz = sizeof (Kde);
  const size_t mcs_sz = sizeof (Mcs) * complexsize.pos;
  const size_t enepara_sz = sizeof (EnePara);
  const size_t temp_sz = sizeof (Temp) * n_tmp;
  const size_t move_scale_sz = sizeof (float) * 6;

  Protein *prt_d[NGPU];
  Psp *psp_d[NGPU];
  Kde *kde_d[NGPU];
  Mcs *mcs_d[NGPU];
  EnePara *enepara_d[NGPU];
  Temp *temp_d[NGPU];
  float *move_scale_d[NGPU];

  for (int i = 0; i < NGPU; ++i) {
    cudaSetDevice (i);
    CUDAMALLOC (prt_d[i], prt_sz, Protein *);
    CUDAMALLOC (psp_d[i], psp_sz, Psp *);
    CUDAMALLOC (kde_d[i], kde_sz, Kde *);
    CUDAMALLOC (mcs_d[i], mcs_sz, Mcs *);
    CUDAMALLOC (enepara_d[i], enepara_sz, EnePara *);
    CUDAMALLOC (temp_d[i], temp_sz, Temp *);
    CUDAMALLOC (move_scale_d[i], move_scale_sz, float *);

    CUDAMEMCPYTOSYMBOL (prt_dc, &prt_d[i], Protein *);
    CUDAMEMCPYTOSYMBOL (psp_dc, &psp_d[i], Psp *);
    CUDAMEMCPYTOSYMBOL (kde_dc, &kde_d[i], Kde *);
    CUDAMEMCPYTOSYMBOL (mcs_dc, &mcs_d[i], Mcs *);
    CUDAMEMCPYTOSYMBOL (enepara_dc, &enepara_d[i], EnePara *);
    CUDAMEMCPYTOSYMBOL (temp_dc, &temp_d[i], Temp *);
    CUDAMEMCPYTOSYMBOL (move_scale_dc, &move_scale_d[i], float *);

    CUDAMEMCPY (prt_d[i], prt, prt_sz, cudaMemcpyHostToDevice);
    CUDAMEMCPY (psp_d[i], psp, psp_sz, cudaMemcpyHostToDevice);
    CUDAMEMCPY (kde_d[i], kde, kde_sz, cudaMemcpyHostToDevice);
    CUDAMEMCPY (mcs_d[i], mcs, mcs_sz, cudaMemcpyHostToDevice);
    CUDAMEMCPY (enepara_d[i], enepara, enepara_sz, cudaMemcpyHostToDevice);
    CUDAMEMCPY (temp_d[i], temp, temp_sz, cudaMemcpyHostToDevice);
    CUDAMEMCPY (move_scale_d[i], &mcpara->move_scale, move_scale_sz, cudaMemcpyHostToDevice);
  }




  // GPU writable arrays that duplicate on multiple GPUs
  const size_t lig_sz = sizeof (Ligand) * n_rep;
  const size_t replica_sz = sizeof (Replica) * n_rep;
  const size_t etotal_sz = sizeof (float) * n_rep;
  const size_t ligmovevector_sz = sizeof (LigMoveVector) * n_rep;
  const size_t acs_temp_exchg_sz = sizeof (int) * n_rep; // acceptance counter
  const size_t ref_matrix_sz = sizeof (ConfusionMatrix);
  //size_t etotal_sz_per_gpu[NGPU];
  //for (int i = 0; i < NGPU; ++i)
  //etotal_sz_per_gpu[i] = sizeof (float) * n_rep_per_gpu[i];

  Ligand *lig_d[NGPU];
  Replica *replica_d[NGPU];
  float *etotal_d[NGPU];
  LigMoveVector *ligmovevector_d[NGPU];
  int *acs_temp_exchg, *acs_temp_exchg_d[NGPU];
  float *ref_matrix_d[NGPU];

  acs_temp_exchg = (int *) malloc (acs_temp_exchg_sz);

  for (int i = 0; i < NGPU; ++i) {
    cudaSetDevice (i);
    CUDAMALLOC (lig_d[i], lig_sz, Ligand *);
    CUDAMALLOC (replica_d[i], replica_sz, Replica *);
    CUDAMALLOC (etotal_d[i], etotal_sz, float *);
    CUDAMALLOC (ligmovevector_d[i], ligmovevector_sz, LigMoveVector *);
    CUDAMALLOC (acs_temp_exchg_d[i], acs_temp_exchg_sz, int *);
    CUDAMALLOC (ref_matrix_d[i], ref_matrix_sz, ConfusionMatrix *);

    CUDAMEMCPYTOSYMBOL (lig_dc, &lig_d[i], Ligand *);
    CUDAMEMCPYTOSYMBOL (replica_dc, &replica_d[i], Replica *);
    CUDAMEMCPYTOSYMBOL (etotal_dc, &etotal_d[i], float *);
    CUDAMEMCPYTOSYMBOL (ligmovevector_dc, &ligmovevector_d[i], LigMoveVector *);
    CUDAMEMCPYTOSYMBOL (acs_temp_exchg_dc, &acs_temp_exchg_d[i], int *);
    CUDAMEMCPYTOSYMBOL (ref_matrix_dc, &ref_matrix_d[i], ConfusionMatrix *);

    CUDAMEMCPY (lig_d[i], lig, lig_sz, cudaMemcpyHostToDevice);
    CUDAMEMCPY (replica_d[i], replica, replica_sz, cudaMemcpyHostToDevice);
  }




  // GPU writable arrays that spreads over multiple GPUs

  // ligrecord[n_rep]
  LigRecord *ligrecord, *ligrecord_d[NGPU];
  const size_t ligrecord_sz = sizeof (LigRecord) * n_rep;
  size_t ligrecord_sz_per_gpu[NGPU];
  for (int i = 0; i < NGPU; ++i)
    ligrecord_sz_per_gpu[i] = sizeof (LigRecord) * n_rep_per_gpu[i];

  ligrecord = (LigRecord *) malloc (ligrecord_sz);
  for (int i = 0; i < NGPU; ++i) {
    cudaSetDevice (i);
    CUDAMALLOC (ligrecord_d[i], ligrecord_sz_per_gpu[i], LigRecord *);
    CUDAMEMCPYTOSYMBOL (ligrecord_dc, &ligrecord_d[i], LigRecord *);
  }

  for (int i = 0; i < NGPU; ++i) {
    cudaGetLastError ();
  }


  // ligand conformation records
  vector < vector < LigRecordSingleStep > > multi_reps_records;

  // launch GPU kernels
  printf ("Start launching kernels\n");

#if NGPU > 1
#include "kernel_cuda_lancher_multigpu.cu"
#elif NGPU == 1
#include "kernel_cuda_lancher_siglegpu.cu"
#endif

  putchar ('\n');
  printf ("================================================================================\n");
  printf ("initial energy state\n");
  printf ("================================================================================\n");
  printf ("rep step vdw ele pmf psp hdb hpc kde lhm dst total\n");
  printf ("0 0");
  LigRecordSingleStep step = multi_reps_records[0][0];
  for (int i = 0; i < MAXWEI; i++)
    printf(" %.3f", step.energy.e[i]);
  putchar ('\n');

  int total_results = multi_reps_records.size();
  SingleRepResult * results = new SingleRepResult[total_results];

  if (!(strlen(mcpara->csv_path) == 0)) {
    printHeader(mcpara);
    vector < vector < LigRecordSingleStep > > :: iterator itr;
    for (itr = multi_reps_records.begin(); itr != multi_reps_records.end(); itr++)
      printStates((*itr), mcpara);
  }
  
  processOneReplica(multi_reps_records[0], &results[0]);

  
  SingleRepResult * first_rep = &results[0];
  printf("================================================================================\n");
  printf("Docking result\n");
  printf("================================================================================\n");
  printf("acceptance ratio\t\t%.3f\n", first_rep->accpt_ratio);
  printf("initial cms\t\t\t%.3f\n", first_rep->init_cms);
  printf("initial rmsd\t\t\t%.3f\n", first_rep->init_rmsd);
  printf("best scored cms\t\t\t%.3f\n", first_rep->best_scored_cms);
  printf("best scored rmsd\t\t%.3f\n", first_rep->best_scored_rmsd);
  printf("best rmsd achieved\t\t%f\n", first_rep->best_achieved_rmsd);
  printf("best cms achieved\t\t%f\n", first_rep->best_achieved_cms);
  printf("pearson between score and rmsd\t%f\n", first_rep->ener_rmsd_p);
  printf("pearson between score and cms\t%f\n", first_rep->ener_cms_p);

  

  delete[]results;

  mclog->ac_mc = 0;
  vector < vector < LigRecordSingleStep > > :: iterator it_rec;
  for (it_rec = multi_reps_records.begin(); it_rec != multi_reps_records.end(); it_rec++)
    mclog->ac_mc += (*it_rec).size();

  // calcuate acceptance counters
  for (int i = 0; i < NGPU; ++i) {
    cudaSetDevice (i);
    CUDAMEMCPY (acs_temp_exchg, acs_temp_exchg_d[i], acs_temp_exchg_sz, cudaMemcpyDeviceToHost);
    for (int j = 0; j < n_rep; ++j) {
      mclog->acs_temp_exchg[j] += acs_temp_exchg[j];
    }
  }
  for (int j = 0; j < n_rep; ++j) {
    mclog->ac_temp_exchg += mclog->acs_temp_exchg[j];
  }








  // free memories
  for (int i = 0; i < NGPU; ++i) {
    cudaSetDevice (i);

    CUDAFREE (curandstate_d[i]);

    CUDAFREE (prt_d[i]);
    CUDAFREE (psp_d[i]);
    CUDAFREE (kde_d[i]);
    CUDAFREE (mcs_d[i]);
    CUDAFREE (enepara_d[i]);
    CUDAFREE (temp_d[i]);
    CUDAFREE (move_scale_d[i]);

    CUDAFREE (lig_d[i]);
    CUDAFREE (replica_d[i]);
    CUDAFREE (etotal_d[i]);
    CUDAFREE (ligmovevector_d[i]);
    CUDAFREE (acs_temp_exchg_d[i]);
    CUDAFREE (ref_matrix_d[i]);

    CUDAFREE (ligrecord_d[i]);
  }

  free (acs_temp_exchg);
  free (ligrecord);
}

