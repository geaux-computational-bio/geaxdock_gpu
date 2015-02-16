#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include "dock.h"
#include "size.h"
#include "load.h"
#include "data.h"
#include "debug.h"


using namespace std;


int
main (int argc, char **argv)
{

#if 0
  if (argc != 3) {
    printf ("usaage: %s [lig/prt] inputfile\n", argv[0]);
    exit (999);
  }

  char mode[8];
  char input[200];
  strcpy (mode, argv[1]);
  strcpy (input, argv[2]);
#endif

  /* load ligand file */
  LigandFile *lig_file = new LigandFile[1];
  ProteinFile *prt_file = new ProteinFile[1];
  Ligand0 *lig = new Ligand0[MAXEN2];
  Protein0 *prt = new Protein0[MAXEN1];

  /* pass args */
  std::string mode = argv[1];

  /* load protein file */
  if (mode == "prt") {
    prt_file->path = argv[2];
    loadProtein (prt_file, prt);
    // cout << "id," << "conf#," << "pts#," << "residue#" << endl;
    cout << prt_file->path << "," << prt_file->conf_total
      << "," << prt_file->pnp << "," << prt_file->pnr << endl;
  }
  /* load ligand file */
  else if (mode == "lig") {
    lig_file->path = argv[2];
    lig_file->molid = "MOLID";
    // loadLigConf ( lig_file );
    loadLigand (lig_file, lig);
    // cout << "id," << "conf#" << "pts#"<< endl;
    cout << lig_file->id << "," << lig_file->raw_conf << "," << lig_file->
      conf_total << ", " << lig_file->lna << endl;
  }
  else
    cout << "give me sth to read!!" << endl;
/*
  for ( int i = 0 ; i < lig->lna; i++)
  {
		  DEBUG_1_("type: ", lig->t[i]);
		  DEBUG_1_("charge: ", lig->c[i]);
  }
*/


  delete[]prt;
  delete[]prt_file;
  delete[]lig;
  delete[]lig_file;

#if 0
  if (strcmp (mode, "prt") == 0) {
    Protein0 *prt = new Protein0[MAXEN1];
    loadProtein (input, prt);
    printf ("%s\t%s\t%d\t%d\n", argv[2], argv[1], prt->pnp, prt->pnr);

    delete[]prt;
  }
#endif


  return 0;
}
