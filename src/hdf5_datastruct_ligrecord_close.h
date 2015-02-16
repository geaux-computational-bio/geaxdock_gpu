//#ifndef HDF5_DATASTRCUT_LIGRECORD_CLOSE_H
//#define HDF5D_ATASTRCUT_LIGRECORD_CLOSE_H


//status = H5Tclose (sz_MAXLIG_t);
  status = H5Tclose (sz_3_t);
  status = H5Tclose (sz_6_t);
  status = H5Tclose (sz_MAXWEI_t);
  status = H5Tclose (sz_LigRecord_t);
  
  status = H5Tclose (Replica_t);
//status = H5Tclose (LigCoord_t);
  status = H5Tclose (Energy_t);
  status = H5Tclose (LigRecordSingleStep_t);
  status = H5Tclose (LigRecord_t);


//#endif
