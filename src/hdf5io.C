#include <cstdlib>
#include <cstdio>
#include <hdf5.h>

#include "dock.h"
#include "size.h"
#include "hdf5io.h"


void
DumpLigRecord (const LigRecord * ligrecord, const int recsize, const char *h5file)
{
  herr_t status;

#include "hdf5_datastruct_ligrecord_create.h"

  // file I/O
  hid_t file = H5Fcreate (h5file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  hid_t space = H5Screate_simple (1, sz_recsize, NULL);
  hid_t dset = H5Dcreate (file, DATASET, LigRecord_t, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite (dset, LigRecord_t, H5S_ALL, H5S_ALL, H5P_DEFAULT, ligrecord);

  // clean data structures
  status = H5Dclose (dset);
  status = H5Sclose (space);
  status = H5Fclose (file);

#include "hdf5_datastruct_ligrecord_close.h"
}






void
ReadLigRecord (LigRecord * ligrecord, const int recsize, const char *h5file)
{
  herr_t status;

#include "hdf5_datastruct_ligrecord_create.h"

  hid_t file = H5Fopen (h5file, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0) {
    printf ("Data file %s do not exist!\n", h5file);
    exit (0);
  }
  hid_t dset = H5Dopen2 (file, "/dset", H5P_DEFAULT);
  status = H5Dread (dset, LigRecord_t, H5S_ALL, H5S_ALL, H5P_DEFAULT, ligrecord);

  // clean data structures
  status = H5Dclose (dset);
  status = H5Fclose (file);

#include "hdf5_datastruct_ligrecord_close.h"
}


