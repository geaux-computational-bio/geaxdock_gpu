#include <stdio.h>
#include <string>
#include <vector>

#include "size.h"
#include "dock.h"
#include "load.h"
#include "util.h"

#include "gtest/gtest.h"
#include "gtest/internal/gtest-internal.h"

using namespace std;

TEST (original_load, 1a07C1)
{
  InputFiles *inputfiles = new InputFiles[1];
  inputfiles->lig_file.id = "1a07C1";
  inputfiles->lig_file.path = "../data/1a07C1/1a07C1.sdf";
  inputfiles->lig_file.molid = "MOLID";
  inputfiles->lhm_file.path = "../data/1a07C1/1a07C1.ff";

  Ligand0 *lig0 = new Ligand0[MAXEN2];
  loadLigand_bk (&inputfiles->lig_file, lig0);

  Ligand0 *lig1 = new Ligand0[MAXEN2];
  loadLigand(&inputfiles->lig_file, lig1);
  trimLigand(inputfiles, lig1);
  
  /*
  for (int i = 0; i < 3; i++) {
    float left = lig0->coord_orig.center[i];
    float right = lig1->coord_orig.center[i];
    cout << left << " " << right << endl;
    const testing::internal::FloatingPoint<float> lhs(left), rhs(right);
    EXPECT_TRUE(lhs.AlmostEquals(rhs));
  }
  
  for (int i = 0; i < 3; i++) {
    float left = lig0->pocket_center[i];
    float right = lig1->pocket_center[i];
    cout << left << " " << right << endl;
    const testing::internal::FloatingPoint<float> lhs(left), rhs(right);
    EXPECT_TRUE(lhs.AlmostEquals(rhs));
  }
  */


  delete[]lig1;
  delete[]lig0;
  delete[]inputfiles;
}


TEST (Load_Ligands, 1a07C1)
{
  string sdf_path = "../data/1a07C1/1a07C1.sdf";
  vector < vector < string > > sections = readLigandSections(sdf_path);

  // should have 1 ligands
  EXPECT_EQ(1, sections.size());

  // load number of ensembles
  vector < string > sect_0 = sections.at(0);
  int ens_total = getLigEnsembleTotal(sect_0);
  // EXPECT_EQ(50, ens_total);

  vector < string > coords;
  coords = getLigEnsembleCoords(sect_0);
  EXPECT_EQ(ens_total, coords.size());

  vector < float > rmsds = getLigEnsembleRmsd(sect_0);
  EXPECT_EQ(ens_total, rmsds.size());

  Ligand0 *lig0 = new Ligand0[MAXEN2];

  int tot_conf = loadOneLigand(sections.at(0), lig0);
  EXPECT_EQ(21, tot_conf);

  float left = 45.6740;
  float right = lig0->coord_orig.x[0];
  const testing::internal::FloatingPoint<float> lhs(left), rhs(right);
  EXPECT_TRUE(lhs.AlmostEquals(rhs));

  Ligand0 * lig1 = &lig0[2];
  left = -5.2246;
  right = lig1->coord_orig.y[0];
  const testing::internal::FloatingPoint<float> lhs1(left), rhs1(right);
  EXPECT_TRUE(lhs1.AlmostEquals(rhs1));

  

  delete[]lig0;

}

TEST (Load_Ligands, 1b9vA)
{
  string sdf_path = "../data/edud/1b9v_4.sdf";
  vector < vector < string > > sections = readLigandSections(sdf_path);

  // should have 6 ligands
  EXPECT_EQ(6, sections.size());

  // load number of ensembles
  vector < string > sect_0 = sections.at(0);
  int ens_total = getLigEnsembleTotal(sect_0);
  EXPECT_EQ(50, ens_total);

  vector < string > coords;
  coords = getLigEnsembleCoords(sect_0);
  EXPECT_EQ(ens_total, coords.size());

  vector < float > rmsds = getLigEnsembleRmsd(sect_0);
  EXPECT_EQ(ens_total, rmsds.size());

  Ligand0 *lig0 = new Ligand0[MAXEN2];

  writeLigProperty (sect_0, lig0);
  EXPECT_EQ(23, lig0->lna);
  EXPECT_EQ(23, lig0->lnb);

  // TODO ATOMIC_TYPES are missing in the edud sdf file TODO
  writeLigAtomProperty (sect_0, lig0);
  
  writeDefaultLigAtomCoord (sect_0, lig0);
  float left = -0.0184;
  float right = lig0->coord_orig.x[1];
  const testing::internal::FloatingPoint<float> lhs(left), rhs(right);
  EXPECT_TRUE(lhs.AlmostEquals(rhs));

  Ligand0 * lig1 = &lig0[1];
  writeEnsAtomCoord (coords.at(0), lig1);
  left = 2.1011;
  right = lig1->coord_orig.x[1];
  const testing::internal::FloatingPoint<float> lhs1(left), rhs1(right);
  EXPECT_TRUE(lhs1.AlmostEquals(rhs1));

  writeEnsAtomCoord (coords.at(2), lig1);
  left = -2.2202;
  right = lig1->coord_orig.z[1];
  const testing::internal::FloatingPoint<float> lhs2(left), rhs2(right);
  EXPECT_TRUE(lhs2.AlmostEquals(rhs2));

  
  moveLigand2ItsCenterFrame(lig0);
  left = -3.6306739;
  right = lig0->coord_orig.x[0];
  const testing::internal::FloatingPoint<float> lhs3(left), rhs3(right);
  EXPECT_TRUE(lhs3.AlmostEquals(rhs3));

  float pocket_center[3];
  loadPocketCenter("../data/edud/1b9v_4.ff", pocket_center);

  left = -9.590;
  right = pocket_center[1];
  const testing::internal::FloatingPoint<float> lhs4(left), rhs4(right);
  EXPECT_TRUE(lhs4.AlmostEquals(rhs4));

  moveLigand2PocketCenter(pocket_center, lig0);

  delete[]lig0;

}


TEST (Load_Ligands, 1e66)
{
  string sdf_path = "../data/edud/1e66_4.sdf";
  vector < vector < string > > sections = readLigandSections(sdf_path);

  // should have 6 ligands
  EXPECT_EQ(6, sections.size());

  // load number of ensembles
  vector < string > sect_0 = sections.at(0);
  int ens_total = getLigEnsembleTotal(sect_0);
  EXPECT_EQ(50, ens_total);

  vector < string > coords;
  coords = getLigEnsembleCoords(sect_0);
  EXPECT_EQ(ens_total, coords.size());

  vector < float > rmsds = getLigEnsembleRmsd(sect_0);
  EXPECT_EQ(ens_total, rmsds.size());
  EXPECT_EQ(ens_total, rmsds.size());

  Ligand0 *lig0 = new Ligand0[MAXEN2];

  writeLigProperty(sect_0, lig0);
  EXPECT_EQ(28, lig0->lna);
  EXPECT_EQ(31, lig0->lnb);


  delete[]lig0;

}
