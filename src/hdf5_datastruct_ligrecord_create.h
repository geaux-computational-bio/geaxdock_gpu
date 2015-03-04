//#ifndef HDF5_DATASTRCUT_LIGRECORD_CREATE_H
//#define HDF5D_ATASTRCUT_LIGRECORD_CREATE_H


/*
struct LigRecord
{
  LigRecordSingleStep step[STEPS_PER_DUMP];
};

struct LigRecordSingleStep
{
  Replica replica;
  Energy energy;
  float movematrix[6]; // // translation x y z, rotation x y z
  int step;
};

struct Replica
{
  int idx_rep; // n_rep, replica
  int idx_tmp; // n_tmp, temperature
  int idx_lig; // n_lig, ligand
  int idx_prt; // n_prt, protein
};

struct Energy
{
  float e[MAXWEI];
  float cms;
  float rmsd;
};
*/



//hsize_t sz_MAXLIG[1] = { MAXLIG };
  hsize_t sz_1[1] = { 1 };
  hsize_t sz_3[1] = { 3 };
  hsize_t sz_6[1] = { 6 };
  hsize_t sz_MAXWEI[1] = { MAXWEI };
  hsize_t sz_STEPS_PER_DUMP[1] = { STEPS_PER_DUMP };
  hsize_t sz_recsize[1] = { recsize };

//hid_t sz_MAXLIG_t = H5Tarray_create (H5T_NATIVE_FLOAT, 1, sz_MAXLIG);
  hid_t sz_1_t = H5Tarray_create (H5T_NATIVE_FLOAT, 1, sz_1);
  hid_t sz_3_t = H5Tarray_create (H5T_NATIVE_FLOAT, 1, sz_3);
  hid_t sz_6_t = H5Tarray_create (H5T_NATIVE_FLOAT, 1, sz_6);
  hid_t sz_MAXWEI_t = H5Tarray_create (H5T_NATIVE_FLOAT, 1, sz_MAXWEI);


  // declare data structures
  hid_t Replica_t = H5Tcreate (H5T_COMPOUND, sizeof (Replica));
  status = H5Tinsert (Replica_t, "idx_rep", HOFFSET (Replica, idx_rep), H5T_NATIVE_INT);
  status = H5Tinsert (Replica_t, "idx_tmp", HOFFSET (Replica, idx_tmp), H5T_NATIVE_INT);
  status = H5Tinsert (Replica_t, "idx_lig", HOFFSET (Replica, idx_lig), H5T_NATIVE_INT);
  status = H5Tinsert (Replica_t, "idx_prt", HOFFSET (Replica, idx_prt), H5T_NATIVE_INT);

/*
  hid_t LigCoord_t = H5Tcreate (H5T_COMPOUND, sizeof (LigCoord));
  status = H5Tinsert (LigCoord_t, "x", HOFFSET (LigCoord, x), sz_MAXLIG_t);
  status = H5Tinsert (LigCoord_t, "y", HOFFSET (LigCoord, y), sz_MAXLIG_t);
  status = H5Tinsert (LigCoord_t, "z", HOFFSET (LigCoord, z), sz_MAXLIG_t);
  status = H5Tinsert (LigCoord_t, "center", HOFFSET (LigCoord, center), sz_3_t);
*/

  hid_t Energy_t = H5Tcreate (H5T_COMPOUND, sizeof (Energy));
  status = H5Tinsert (Energy_t, "e", HOFFSET (Energy, e), sz_MAXWEI_t);
  status = H5Tinsert (Energy_t, "cms", HOFFSET (Energy, cms), sz_1_t);
  status = H5Tinsert (Energy_t, "rmsd", HOFFSET (Energy, rmsd), sz_1_t);

  hid_t LigRecordSingleStep_t = H5Tcreate (H5T_COMPOUND, sizeof (LigRecordSingleStep));
  status = H5Tinsert (LigRecordSingleStep_t, "replica", HOFFSET (LigRecordSingleStep, replica), Replica_t);
  status = H5Tinsert (LigRecordSingleStep_t, "energy", HOFFSET (LigRecordSingleStep, energy), Energy_t);
//status = H5Tinsert (LigRecordSingleStep_t, "coord", HOFFSET (LigRecordSingleStep, coord), LigCoord_t);
  status = H5Tinsert (LigRecordSingleStep_t, "movematrix", HOFFSET (LigRecordSingleStep, movematrix), sz_6_t);
  status = H5Tinsert (LigRecordSingleStep_t, "step", HOFFSET (LigRecordSingleStep, step), H5T_NATIVE_INT);




  hid_t sz_LigRecord_t = H5Tarray_create (LigRecordSingleStep_t, 1, sz_STEPS_PER_DUMP);
  hid_t LigRecord_t = H5Tcreate (H5T_COMPOUND, sizeof (LigRecord));
  status = H5Tinsert (LigRecord_t, "step", HOFFSET (LigRecord, step), sz_LigRecord_t);


//#endif

