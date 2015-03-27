#include <iostream>
#include "kgs.h"

using namespace std;

map < int, vector < int > >
GetClusters(int* clusterid, int ncluster, int nobj)
{
  map < int, vector < int > > clusters;
  for (int i = 0; i < ncluster; i++)
    for (int j = 0; j < nobj; j++)
      {
        if (clusterid[j] == i)
          clusters[i].push_back(j);
      }
  return clusters;
}

map < int, double >
Distances2Others(vector < int > & members, double** distmatrix)
{
  vector < int > :: iterator itm1;
  vector < int > :: iterator itm2;

  map < int, double > dists;
  for (itm1 = members.begin(); itm1 != members.end(); itm1 ++) {
    int my_idx = (*itm1);
    double tot_dist_between = 0.0;

    for (itm2 = members.begin(); itm2 != members.end(); itm2 ++) {
      int other_idx = (*itm2);
      double dist_between = distmatrix[my_idx][other_idx];
      tot_dist_between += dist_between;
    }

    dists[my_idx] = tot_dist_between;
  }
  return dists;
}

int
FindMedoid(map < int, double > & pt_and_its_dist_to_others)
{

  map < int, double > :: iterator itp;
  int medoid_idx = -1;
  double min_dist = ( (double) pt_and_its_dist_to_others.size() ) * MAX_DIST;

  for (itp = pt_and_its_dist_to_others.begin();
       itp != pt_and_its_dist_to_others.end();
       ++itp)
    {
      int my_idx = itp->first;
      double dist = itp->second;
      if (dist < min_dist) {
        min_dist = dist;
        medoid_idx = my_idx;
      }
    }

  return medoid_idx;
}

double
SpreadOfCluster(map < int, double > &pt_and_its_dist_to_others)
{
  assert(pt_and_its_dist_to_others.size() > 1);

  map < int, double > :: iterator itp;
  double sum = 0.0;
  
  for (itp = pt_and_its_dist_to_others.begin();
       itp != pt_and_its_dist_to_others.end();
       ++itp)
    {
      double dist = itp->second;
      sum += dist;
    }

  int sz = pt_and_its_dist_to_others.size();
  int counts = sz * (sz - 1);
  double spread = sum / (double) counts;
  return spread;
}

double
AveSpread(map < int, vector < int > > & clusters, double** dist_matrix)
{
  map < int, map < int, double > > dist_to_others;
  map < int , vector < int > > :: iterator itc;
  int non_outliers = 0;
  double tot_spread = 0.0;
  for (itc = clusters.begin(); itc != clusters.end(); itc ++)
    {
    map < int, double > pt_and_its_dist_to_others  = Distances2Others(itc->second, dist_matrix);
    if (pt_and_its_dist_to_others.size() > 1)
      {
        double spread = SpreadOfCluster(pt_and_its_dist_to_others);
        tot_spread += spread;
        non_outliers += 1;
      }
    }
  assert (non_outliers != 0);
  return tot_spread / (double) non_outliers;
}

int
KGS(Node* tree, int* clusterid, double** dist_matrix, int nobj, bool show_penalties)
{
  int ncluster;
  int max_clusters = nobj - 1;
  vector < double > ave_spreads;
  map < int , vector < int > > clusters;

  for (ncluster = max_clusters; ncluster != 1; --ncluster)
    {
      cuttree (nobj, tree, ncluster, clusterid);
      clusters = GetClusters(clusterid, ncluster, nobj);
      double ave_spread = AveSpread(clusters, dist_matrix);
      ave_spreads.push_back(ave_spread);
      // cout << "# clusters " << ncluster << " average spread " << ave_spread << endl;
    }

  double max_ave_spread = *std::max_element(ave_spreads.begin(), ave_spreads.end());
  double min_ave_spread = *std::min_element(ave_spreads.begin(), ave_spreads.end());
  double multiplier = (nobj - 2) / (max_ave_spread - min_ave_spread);

  // printf("--------------------------------------------------------------------------------\n");
  vector < double > :: iterator its;
  vector < double > penalties;
  for (its = ave_spreads.begin(), ncluster = max_clusters;
       its != ave_spreads.end();
       ++its, --ncluster)
    {
      double my_spread = *its;
      double normed_spread = multiplier * (my_spread - min_ave_spread) + 1;
      double penalty = normed_spread + ncluster;
      penalties.push_back(penalty);
      if (show_penalties)
        cout << "# clusters " << ncluster << " penalty " << penalty << endl;
    }

  // printf("--------------------------------------------------------------------------------\n");
  its = penalties.begin();
  int lowest_penalty_cluster_num = nobj - 1;
  double lowest_penalty = (*its);
  for (its = penalties.begin(), ncluster = max_clusters;
       its != penalties.end();
       ++its, --ncluster)
    {
      if (*its < lowest_penalty)
        {
          lowest_penalty = *its;
          lowest_penalty_cluster_num = ncluster;
        }
    }

  // cout << "with cluster # " << lowest_penalty_cluster_num << " has lowest penalty value: " << lowest_penalty << endl;

  // printf("--------------------------------------------------------------------------------\n");
  return lowest_penalty_cluster_num;
  
}
