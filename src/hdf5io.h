#ifndef  HDF5IO_H
#define  HDF5IO_H


#define DATASET "dset"

void
DumpLigRecord (const LigRecord * ligrecord, const int recsize, const char *h5file);

void
ReadLigRecord (LigRecord * ligrecord, const int recsize, const char *h5file);



#endif
