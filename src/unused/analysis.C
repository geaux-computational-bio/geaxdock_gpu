#include <cstdlib>
#include <cstdio>

#include "dock.h"
#include "util.h"
#include "size.h"
#include "hdf5io.h"


int
main (int argc, char **argv)
{
  if (argc < 2) {
    fprintf (stderr, "usage: %s  <input file>\n", argv[0]);
    printf ("-nl <number of show lines>\n");
    printf ("-l <ligand conf number>\n");
    printf ("-p <protein conf number>\n");
  }

  // default settings
  int num_show_line = 0;
  int lig_conf_num = 0;
  int prt_conf_num = 0;
  int show_energy = 0;
  int show_rep = 0;
  int myreplica = 0;

  for ( int i = 0; i < argc; i++ ) {
    if ( !strcmp(argv[i],"-nl")  && i < argc ) 
      num_show_line = atoi(argv[i+1]);
    if ( !strcmp(argv[i],"-l")  && i < argc ) 
      lig_conf_num = atoi(argv[i+1]);
    if ( !strcmp(argv[i],"-e")  && i < argc ) 
      show_energy = 1;
    if ( !strcmp(argv[i],"-r")  && i < argc ) 
      show_rep = 1;
    if ( !strcmp(argv[i],"-p")  && i < argc ) 
      prt_conf_num = atoi(argv[i+1]);
    if ( !strcmp(argv[i],"-rep")  && i < argc ) 
      myreplica = atoi(argv[i+1]);
  }

  ComplexSize complexsize;
  complexsize.n_prt = 3;
  complexsize.n_tmp = MAXTMP;
  complexsize.n_lig = 20;
  complexsize.n_rep = complexsize.n_lig * complexsize.n_prt * complexsize.n_tmp;
  complexsize.n_pos = 0; // unused, the value does not matter

  LigRecord *ligrecord;
  size_t ligrecord_sz = sizeof (LigRecord) * complexsize.n_rep;
  ligrecord = (LigRecord *) malloc (ligrecord_sz);
  ReadLigRecord (ligrecord, complexsize.n_rep, argv[argc-1]);

  // int repp_begin = 0;
  // int repp_end = 22;
  int iter_begin = 0;
  // int iter_end = STEPS_PER_DUMP - 1;
  int iter_end = minimal_int (STEPS_PER_DUMP, num_show_line) - 1;
  int arg = 2;

  if (show_energy == 1)
    PrintLigRecord (ligrecord, STEPS_PER_DUMP, myreplica, iter_begin, iter_end, arg);
  //PrintRepRecord (ligrecord, STEPS_PER_DUMP, repp_begin, repp_end, iter_begin, iter_end, arg);
  if (show_rep == 1)
    PrintRepRecord2 (ligrecord, complexsize, STEPS_PER_DUMP, 
		     lig_conf_num, prt_conf_num, 
		     iter_begin, iter_end, arg);
  //PrintMoveRecord (ligrecord, STEPS_PER_DUMP, myreplica, iter_begin, iter_end, arg);

  free (ligrecord);
  return 0;
}

