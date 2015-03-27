#ifndef KGS_H
#define KGS_H

#include <map>
#include <vector>
#include <algorithm>

#include <assert.h>

#include "size.h"

extern "C" {
#include "./modules/cluster-1.52a/src/cluster.h" /* The C Clustering Library */
}

using namespace std;

// map between each cluster to its members
map < int, vector < int > > GetClusters(int* clusterid, int ncluster, int nobj);

// dictionary of a point and its total distance to others
map < int, double > Distances2Others(vector < int > &members, double** distmatrix);

// find the idx of the medoid in a cluster
// defined as the point with the minimum average distance to all others
int FindMedoid(map < int, double > & pt_and_its_dist_to_others);

double SpreadOfCluster(map < int, double > &pt_and_its_dist_to_others);

double AveSpread(map < int, vector < int > > & clusters, double** dist_matrix);

// find the number of clusters that gives the smallest KGS penalty
int KGS(Node* tree, int* clusterid, double** dist_matrix, int nobj, bool show_penalties);

#endif // UTIL_H
