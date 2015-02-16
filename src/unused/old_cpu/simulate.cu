// ...

#if 0
  MonteCarlo_Init (lig, prt, psp, kde, mcs, enepara);

  // PrintLigCoord (lig, 0);
  // PrintLigCoord (lig, 1);

  double  t0 = host_time_now ();
  MonteCarlo (lig, prt, psp, kde, mcs, enepara, mcpara);
  double  t1 = host_time_now ();
  printf ("ms per MC sweep:\t%7.3f\n", (t1 - t0) * 1000000 / mcpara->step);
#endif

// ...
