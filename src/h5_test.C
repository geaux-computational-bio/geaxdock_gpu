#include <stdio.h>
#include <string>
#include <vector>

#include "size.h"
#include "dock.h"
#include "load.h"
#include "util.h"
#include "hdf5io.h"
#include "hdf5io.h"

#include "gtest/gtest.h"
#include "gtest/internal/gtest-internal.h"

using namespace std;

TEST (load_h5, 1a07C1)
{
  InputFiles *inputfiles = new InputFiles[1];
  Ligand0 *lig0 = new Ligand0[MAXEN2];
  Protein0 *prt0 = new Protein0[MAXEN1];

  inputfiles->lig_file.id = "1a07C1";
  inputfiles->lig_file.path = "../data/1a07C1/1a07C1.sdf";
  inputfiles->lig_file.molid = "MOLID";
  inputfiles->lhm_file.path = "../data/1a07C1/1a07C1-0.8.ff";
  inputfiles->prt_file.path = "../data/1a07C1/1a07C.pdb";

  loadLigand(inputfiles, lig0);
  loadProtein (&inputfiles->prt_file, prt0);

  EXPECT_EQ(3, inputfiles->prt_file.conf_total);

  int n_tmp = 1;
  int n_lig = inputfiles->lig_file.conf_total;
  int n_prt = inputfiles->prt_file.conf_total;
  int n_rep = n_tmp * n_lig * n_prt;

  LigRecord * ligrecord;
  size_t ligrecord_sz = sizeof (LigRecord) * n_rep;
  ligrecord = (LigRecord *) malloc (ligrecord_sz);

  string h5_path = "../data/output_20150224_204013/a_0000.h5";
  const char * path = h5_path.c_str();

  ReadLigRecord(ligrecord, n_rep, path);

  // LigRecordSingleStep step = ligrecord->step[0];
  // for (int i = 0; i < STEPS_PER_DUMP; i++) {
  for (int i = 0; i < 100; i++) {
    cout << ligrecord->step[i].step << endl;
  }


  free (ligrecord);
    

  delete[]lig0;
  delete[]prt0;
  delete[]inputfiles;
}
