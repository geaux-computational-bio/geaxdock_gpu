#ifndef  RUN_H
#define  RUN_H

#include <map>
#include <vector>
#include "dock.h"

using namespace std;

void
Run (const Ligand *,
     const Protein *,
     const Psp *,
     const Kde *,
     const Mcs *,
     const EnePara *,
     const Temp *,
     const Replica *,
     const McPara *,
     McLog *,
     map < int, vector < LigRecordSingleStep > > & multi_reps_records,
     const ComplexSize);


#endif // RUN_H


