



double t1 = HostTimeNow ();

int i = 0;
cudaSetDevice (i);
CUDAKERNELSYNC (ResetCounter_d, dim_grid, dim_block, rep_begin[i], rep_end[i]);
CUDAKERNELSYNC (MonteCarlo_Init_d, dim_grid, dim_block, rep_begin[i], rep_end[i]);

int s1 = 0;
int est_tot_rec = mcpara->steps_per_dump * complexsize.n_rep;
printf("estimated total records: %d\n", est_tot_rec);
// int est_tot_rec = MINIMUM_REC;

while(CountValidRecords(multi_reps_records) < est_tot_rec) {

  // skip the reset if since it has been performed before MC_init
  if (s1 != 0)
    CUDAKERNELSYNC (ResetCounter_d, dim_grid, dim_block, rep_begin[i], rep_end[i]);

  double t0 = HostTimeNow ();

  for (int s2 = 0; s2 < mcpara->steps_per_dump; s2 += mcpara->steps_per_exchange) {
    CUDAKERNELSYNC (MonteCarlo_d, dim_grid, dim_block, rep_begin[i], rep_end[i], s1, s2);
# if IS_EXCHANGE == 1
    const int mode_l = 4; // ligand exchange mode
    const int mode_t = !((s2 / mcpara->steps_per_exchange) % 2); // temperature exchange mode
    CUDAKERNELSYNC (ExchangeReplicas_d, dim_grid, dim_block, mode_l, mode_t);
# endif
  }

  // accumulate for compute time
  mclog->t0 += HostTimeNow () - t0;


  // copy ligand record from GPU to CPU memory
  CUDAMEMCPY (&ligrecord[rep_begin[i]], ligrecord_d[i], ligrecord_sz_per_gpu[i], cudaMemcpyDeviceToHost);

  // gather ligand record
  for (int rep = rep_begin[i]; rep <= rep_end[i]; ++rep) {
    for (int s = 0; s < ligrecord[rep].next_ptr; ++s) {
      LigRecordSingleStep my_step = ligrecord[rep].step[s];
      multi_reps_records[rep].push_back(my_step);
    }
  }

  printf("# points\t\t\t%d\n", CountValidRecords(multi_reps_records));
  fflush (stdout);


#if IS_OUTPUT == 0
    
  // dump ligand record from CPU memory to disk
  if (!(strlen(mcpara->hdf_path) == 0))
    DumpLigRecord (ligrecord, n_rep, mcpara->hdf_path);

#endif
    
  // accumulate for wall time (compute time plus I/O time)
  s1 += mcpara->steps_per_dump;
 }

mclog->t1 += HostTimeNow () - t1;
mclog->steps_total = s1;

int trials = complexsize.n_rep * s1;
mclog->ar = (float) CountValidRecords(multi_reps_records) / (float) trials;


