

// no replica exchange



for (int i = 0; i < NGPU; ++i) {
  cudaSetDevice (i);
  CUDAKERNELASYNC (MonteCarlo_Init_d, dim_grid, dim_block, rep_begin[i], rep_end[i]);
 }

for (int s1 = 0; s1 < mcpara->steps_total; s1 += mcpara->steps_per_dump) {

  double t0 = HostTimeNow ();
  printf ("\t%d / %d \r", s1, mcpara->steps_total);
  fflush (stdout);

  for (int s2 = 0; s2 < mcpara->steps_per_dump; s2 += mcpara->steps_per_exchange) {
    for (int i = 0; i < NGPU; ++i) {
      cudaSetDevice (i);
      CUDAKERNELASYNC (MonteCarlo_d, dim_grid, dim_block, rep_begin[i], rep_end[i], s1, s2);
    }
  }


  cudaDeviceSynchronize ();
  // accumulate for compute time
  mclog->t0 += HostTimeNow () - t0;

#if IS_OUTPUT == 1
  for (int i = 0; i < NGPU; ++i) {
    cudaSetDevice (i);
    CUDAMEMCPY (&ligrecord[rep_begin[i]], ligrecord_d[i], ligrecord_sz_per_gpu[i], cudaMemcpyDeviceToHost);
  }

  // dump ligand record from CPU memory to disk
  char myoutputfile[MAXSTRINGLENG];
  sprintf(myoutputfile, "%s/%s_%04d.h5", mcpara->outputdir, mcpara->outputfile, s1 / mcpara->steps_per_dump);
  DumpLigRecord (ligrecord, n_rep, myoutputfile);
#endif

  mclog->t1 += HostTimeNow () - t0;
}






#if 0

cudaStream_t stream[NGPU];
for (int i = 0; i < NGPU; ++i)
  cudaStreamCreate (&stream[i]);

for (int i = 0; i < NGPU; ++i) {
  cudaSetDevice (i);
  CUDAKERNELSTREAMSYNC (MonteCarlo_Init_d, dim_grid, dim_block, 0, stream[i], rep_begin[i], rep_end[i]);
 }

for (int i = 0; i < NGPU; ++i)
  cudaStreamSynchronize (stream[i]);



for (int s1 = 0; s1 < mcpara->steps_total; s1 += mcpara->steps_per_dump) {

  double t0 = HostTimeNow ();
  printf ("\t%d / %d \r", s1, mcpara->steps_total);
  fflush (stdout);

  for (int s2 = 0; s2 < mcpara->steps_per_dump; s2 += mcpara->steps_per_exchange) {

    // monte carlo
    for (int i = 0; i < NGPU; ++i) {
      cudaSetDevice (i);
      CUDAKERNELSTREAM (MonteCarlo_d, dim_grid, dim_block, 0, stream[i], rep_begin[i], rep_end[i], s1, s2);
    }

    /*
    // exchange
    // gather to GPU0, and then scater from GPU0
    for (int i = 1; i < NGPU; ++i)
      cudaMemcpyPeerAsync(etotal_d[0], 0, etotal_d[i], i, etotal_sz, stream[0]);
  
    // may duplicate computation and eliminate the following transformation
    cudaSetDevice (0);
    ExchangeReplicas_d <<<dim_grid, dim_block, 0, stream[0] >>> ();
    
    for (int i = 1; i < NGPU; ++i)
      cudaMemcpyPeerAsync(replica_d[i], i, replica_d[0], 0, replica_sz, stream[i]);
  
    for (int i = 0; i < NGPU; ++i)
      cudaStreamSynchronize (stream[i]);
    */
  }


  // accumulate for compute time
  mclog->t0 += HostTimeNow () - t0;


  // dump: GPU -> CPU -> disk
  /*
#if IS_OUTPUT == 1
  // copy ligand record from GPU to CPU memory
  for (int i = 0; i < NGPU; ++i) {
    cudaSetDevice (i);
    CUDAMEMCPY (&ligrecord[rep_begin[i]], ligrecord_d[i], ligrecord_sz_per_gpu[i], cudaMemcpyDeviceToHost);
  }

  // dump ligand record from CPU memory to disk
  char myoutputfile[MAXSTRINGLENG];
  sprintf(myoutputfile, "%s/%s_%04d.h5", mcpara->outputdir, mcpara->outputfile, s1 / mcpara->steps_per_dump);
  DumpLigRecord (ligrecord, n_rep, myoutputfile);
#endif
  */


  // accumulate for wall time (compute time plus I/O time)
  mclog->t1 += HostTimeNow () - t0;
}



for (int i = 0; i < NGPU; ++i)
  cudaStreamDestroy (stream[i]);


#endif

