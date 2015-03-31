#ifndef POST_MC_H
#define POST_MC_H

#include <vector>
#include <map>

#include "dock.h"
#include "util.h"

using namespace std;

void post_mc(map < int, vector < LigRecordSingleStep > > & multi_reps_records,
             Ligand* lig, 
             const Protein* const prt, 
             const EnePara* const enepara, 
             const McPara* const mcpara);




#endif /* POST_MC_H */
