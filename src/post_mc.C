#include <stdio.h>

#include "post_mc.h"


void 
post_mc(map < int, vector < LigRecordSingleStep > > & multi_reps_records,
        Ligand* lig, 
        const Protein* const prt, 
        const EnePara* const enepara, 
        const McPara* const mcpara)
{
  putchar ('\n');
  printf ("================================================================================\n");
  printf ("initial energy state\n");
  printf ("================================================================================\n");
  printf ("rep step vdw ele pmf psp hdb hpc kde lhm dst total\n");
  printf ("0 0");
  LigRecordSingleStep step = multi_reps_records[0][0];
  for (int i = 0; i < MAXWEI; i++)
    printf(" %.3f", step.energy.e[i]);
  putchar ('\n');

  int total_results = multi_reps_records.size();
  SingleRepResult * results = new SingleRepResult[total_results];

  /* print traces */
  // if (!(strlen(mcpara->csv_path) == 0)) {
  //   printHeader(mcpara);
  //   map < int, vector < LigRecordSingleStep > > :: iterator itr;
  //   for (itr = multi_reps_records.begin(); itr != multi_reps_records.end(); itr++)
  //     printStates(itr->second, mcpara);
  // }

  
  /* clustering */
  
  // string clustering_method = "k";
  // string clustering_method = "a";
  string clustering_method = "c";
  vector < Medoid > medoids;
  
  for(size_t i = 0; i < multi_reps_records.size(); ++i) {
    medoids = clusterOneRepResults(multi_reps_records[i],
                                   clustering_method, lig, prt, enepara);
    // if (!(strlen(mcpara->csv_path) == 0)) {
    //   printStates(medoids, mcpara);
    // }
  }


  // processOneReplica(multi_reps_records[0], &results[0]);

  delete[]results;
}
