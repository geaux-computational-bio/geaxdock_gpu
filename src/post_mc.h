#ifndef POST_MC_H
#define POST_MC_H

#include <vector>
#include <map>

#include "dock.h"
#include "util.h"

using namespace std;

std::vector<Medoid>
    post_mc(map<int, vector<LigRecordSingleStep> > multi_reps_records,
            Ligand *lig, const Protein *const prt, const EnePara *const enepara,
            const McPara *const mcpara);

vector<Medoid>
cluster_trajectories(map < int, vector < LigRecordSingleStep > > & multi_reps_records,
                     Ligand* lig, int n_lig,
                     const Protein* const prt, 
                     const EnePara* const enepara);

void opt_ff(map<int, vector<LigRecordSingleStep> > &multi_reps_records,
            Ligand *lig, int n_lig, const Protein *const prt,
            const EnePara *const enepara, const McPara *const mcpara);

#endif /* POST_MC_H */
