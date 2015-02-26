


vector < vector < LigRecordSingleStep > > record_multi_rep;


  int i = 0;
  cudaSetDevice (i);
  CUDAKERNELSYNC (MonteCarlo_Init_d, dim_grid, dim_block, rep_begin[i], rep_end[i]);

  for (int s1 = 0; s1 < mcpara->steps_total; s1 += mcpara->steps_per_dump) {

    double t0 = HostTimeNow ();
    printf ("\t%d / %d \r", s1, mcpara->steps_total);
    fflush (stdout);

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
    for (int rep = rep_begin[i]; rep < rep_end[i]; ++rep) {
      vector < LigRecordSingleStep > record_single_rep;
      for (int s = 0; s < ligrecord[rep].next_ptr; ++s) {
	record_single_rep.push_back (ligrecord[rep].step[s]);
      }
      record_multi_rep.push_back (record_single_rep);
    }
    // now "record_multi_rep" holds everything




#if IS_OUTPUT == 1
    /*
    // dump ligand record from CPU memory to disk
    char myoutputfile[MAXSTRINGLENG];
    sprintf(myoutputfile, "%s/%s_%04d.h5", mcpara->outputdir, mcpara->outputfile, s1 / mcpara->steps_per_dump);
    DumpLigRecord (ligrecord, n_rep, myoutputfile);
    */
#endif
    
    // accumulate for wall time (compute time plus I/O time)
    mclog->t1 += HostTimeNow () - t0;
  }

