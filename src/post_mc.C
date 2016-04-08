#include <stdio.h>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <climits>
#include <cassert>


#include "dock.h"
#include "util.h"
#include "post_mc.h"

template <class Key, class Value>
static
unsigned long mapSize(const std::map<Key,Value> &map){
  unsigned long size = sizeof(map);
  for(typename std::map<Key,Value>::const_iterator it = map.begin(); it != map.end(); ++it){
    size += it->first.size();
    size += it->second.size();
  }
  return size;
}

template <class Value>
static
unsigned long vecSize(const std::vector<Value> &vec)
{
  unsigned long size = sizeof(vec);
  for (auto it = vec.begin(); it != vec.end(); ++it) {
    size += sizeof(*it);
  }
  return size;
}


template <class Key, class Value>
static
unsigned long mapSize(const std::map<Key, std::vector<Value>> & map)
{
  unsigned long size = sizeof(map);
  for(auto it = map.begin(); it != map.end(); ++it)
    size += vecSize(it->second);
  return size;
}


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
  if (!(strlen(mcpara->csv_path) == 0)) {
    printHeader(mcpara);
    map < int, vector < LigRecordSingleStep > > :: iterator itr;
    for (itr = multi_reps_records.begin(); itr != multi_reps_records.end(); itr++)
      printStates(itr->second, mcpara);
  }

  
  /* clustering */
  
  // string clustering_method = "k";
  // string clustering_method = "a";
  // string clustering_method = "c";
  // vector < Medoid > medoids;
  
  // for(size_t i = 0; i < multi_reps_records.size(); ++i) {
  //   medoids = clusterOneRepResults(multi_reps_records[i],
  //                                  clustering_method, lig, prt, enepara);
  //   // if (!(strlen(mcpara->csv_path) == 0)) {
  //   //   printStates(medoids, mcpara);
  //   // }
  // }


  // processOneReplica(multi_reps_records[0], &results[0]);

  delete[]results;
}

vector<Medoid>
cluster_trajectories(map < int, vector < LigRecordSingleStep > > & multi_reps_records,
                     Ligand* lig, int n_lig,
                     const Protein* const prt, 
                     const EnePara* const enepara)
{
  vector < LigRecordSingleStep > records;

  for (auto it = multi_reps_records.begin(); it != multi_reps_records.end(); ++it) {

    records.insert(records.end(), 
                   make_move_iterator(it->second.begin()),
                   make_move_iterator(it->second.end()));

    it->second.clear();
  }

  sort(records.begin(), records.end(), cmsLargerThan);
  // std::random_shuffle(records.begin(), records.end());


  // for (auto it = records.begin(); it != records.end(); ++it) {
  //   LigRecordSingleStep* s = &(*it);
  //   cout << getCMS(s) << endl;
  // }
  
  assert(records.size() > MINIMUM_REC);
  size_t num_grp = 20;
  size_t total_samples = 5000;
  int num_cluster_each_grp = (int) (total_samples / num_grp);
  size_t num_samples_each_grp = MINIMUM_REC / num_grp;
  
  cout << "================================================================================" << endl;
  cout << "clustering" << endl;
  cout << "================================================================================" << endl;
  cout << "# samples\t\t\t" << MINIMUM_REC << endl;
  cout << "# groups\t\t\t" << num_grp << endl;
  cout << "# samples each group\t\t" << num_samples_each_grp << endl;
  cout << "# cluster each group\t\t" << num_cluster_each_grp << endl;
  cout << "clustering ratio\t\t" << (float) num_samples_each_grp / num_cluster_each_grp << endl;

  vector < Medoid > all_medoids;
  
  for (size_t grp_idx = 0; grp_idx < num_grp; ++grp_idx) {
    cout << "clustering group\t\t" << grp_idx << endl;
    vector<LigRecordSingleStep> steps(records.begin() + grp_idx * num_samples_each_grp, 
                                      records.begin() + (grp_idx + 1) * num_samples_each_grp);
    assert(steps.size() > 0);
    assert(steps.size() < INT_MAX);

    vector < Medoid > medoids = clusterCmsByAveLinkage(steps, num_cluster_each_grp, n_lig, lig, prt, enepara);

    all_medoids.insert(all_medoids.end(),
                       make_move_iterator(medoids.begin()),
                       make_move_iterator(medoids.end()));
  }

  return all_medoids;
}

