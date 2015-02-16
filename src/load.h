#ifndef  LOAD_H
#define  LOAD_H

#include <vector>

#include "dock.h"

using namespace std;

void loadTrace(TraceFile *, float *);

void loadLigConf (LigandFile *);
void loadLigand (LigandFile *, Ligand0 *);

void loadPrtConf (ProteinFile *, Protein0 *);
void loadProtein (ProteinFile *, Protein0 *);

void loadLHM (LhmFile *, Psp0 *, Kde0 *, Mcs0 *);
void loadEnePara (EneParaFile *, EnePara0 *);

void loadWeight(WeightFile *, EnePara0 *);
void loadNorPara(NorParaFile *, EnePara0 *);  // load the normalization parameter a and b

vector < vector < float > > read2D(TraceFile *);


#endif // LOAD_H

