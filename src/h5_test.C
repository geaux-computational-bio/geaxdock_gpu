#include <stdio.h>
#include <string>
#include <vector>
#include <cstring>
#include <cmath>

#include "size.h"
#include "dock.h"
#include "load.h"
#include "util.h"
#include "hdf5io.h"
#include "hdf5io.h"

#include "gtest/gtest.h"
#include "gtest/internal/gtest-internal.h"

using namespace std;

float getTotalEner(LigRecordSingleStep step)
{
  return step.energy.e[MAXWEI - 1];
}

TEST (Energy, 1a07C1)
{
  // testing loadLigand and loadLigand_bk
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

  LigRecord * old_ligrecord;
  size_t old_ligrecord_sz = sizeof (LigRecord) * n_rep;
  old_ligrecord = (LigRecord *) malloc (old_ligrecord_sz);

  // using new loading function
  string h5_path = "../data/output_20150227_163950/a_0000.h5";
  const char * path = h5_path.c_str();
  
  // using old loading function
  string old_h5_path = "../data/output_20150227_164804/a_0000.h5";
  const char * old_path = old_h5_path.c_str();


  // test if two loading functions affect energy calculation
  for (int idx_rep = 0; idx_rep < n_rep; idx_rep++)
    {
      vector < LigRecordSingleStep > records;
      ReadLigRecord(ligrecord, n_rep, path);
      int tot_records = checkRedundancy(records, idx_rep, ligrecord);

      vector < LigRecordSingleStep > old_records;
      ReadLigRecord(old_ligrecord, n_rep, old_path);
      int old_tot_records = checkRedundancy(old_records, idx_rep, old_ligrecord);

      float new_ener = getTotalEner(records[0]);
      float old_ener = getTotalEner(old_records[0]);
      assert (fabs(new_ener - old_ener) < 0.01);
      // printf ("# new : %f # old : %f\n", new_ener, old_ener);
      // printf ("# new : %d # old : %d\n", tot_records, old_tot_records);
    }

  free (old_ligrecord);
  free (ligrecord);

  delete[]lig0;
  delete[]prt0;
  delete[]inputfiles;
  
}
