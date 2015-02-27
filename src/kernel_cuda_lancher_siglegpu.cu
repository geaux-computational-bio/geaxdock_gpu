




  int i = 0;
  cudaSetDevice (i);
  CUDAKERNELSYNC (ResetCounter_d, dim_grid, dim_block, rep_begin[i], rep_end[i]);
  CUDAKERNELSYNC (MonteCarlo_Init_d, dim_grid, dim_block, rep_begin[i], rep_end[i]);

  for (int s1 = 0; s1 < mcpara->steps_total; s1 += mcpara->steps_per_dump) {

    double t0 = HostTimeNow ();
    fflush (stdout);

    // the reset has been done before MC_init, skip it
    if (s1 != 0)
      CUDAKERNELSYNC (ResetCounter_d, dim_grid, dim_block, rep_begin[i], rep_end[i]);

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

    int total_size = 0;
    for (int rep = rep_begin[i]; rep <= rep_end[i]; ++rep)
      total_size += ligrecord[rep].next_ptr;

    // printf("total reps in the first GPU:\t%d\n", rep_end[i]);
    printf("%d records in %d replica\n", ligrecord[0].next_ptr, 0);
    printf("%d records in the first GPU:\n", total_size);


    // gather ligand record
    for (int rep = rep_begin[i]; rep <= rep_end[i]; ++rep) {
      vector < LigRecordSingleStep > single_rep_records;
      //printf("%d records in %d replica\n", ligrecord[rep].next_ptr, rep);
      for (int s = 0; s < ligrecord[rep].next_ptr; ++s) {
        single_rep_records.push_back (ligrecord[rep].step[s]);
      }
      multi_reps_records.push_back(single_rep_records);
    }



#if IS_OUTPUT == 1
    
    // dump ligand record from CPU memory to disk
    char myoutputfile[MAXSTRINGLENG];
    sprintf(myoutputfile, "%s/%s_%04d.h5", mcpara->outputdir, mcpara->outputfile, s1 / mcpara->steps_per_dump);
    DumpLigRecord (ligrecord, n_rep, myoutputfile);
#endif
    
    // accumulate for wall time (compute time plus I/O time)
    mclog->t1 += HostTimeNow () - t0;
  }

