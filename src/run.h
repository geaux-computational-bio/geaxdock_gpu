#ifndef  RUN_H
#define  RUN_H

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
     vector < vector < LigRecordSingleStep > > &,
     const ComplexSize);


#endif // RUN_H


