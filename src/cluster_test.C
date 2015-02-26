#include <stdio.h>
#include <string>
#include <vector>
#include <cstring>

#include "size.h"
#include "dock.h"
#include "load.h"
#include "util.h"
#include "hdf5io.h"
#include "hdf5io.h"

#include "gtest/gtest.h"
#include "gtest/internal/gtest-internal.h"

extern "C" {
#include "kmeans.h"
}

using namespace std;


TEST(Kmeans, run_sample)
{
  int     isBinaryFile, is_output_timing;
  
  int     numClusters, numCoords, numObjs;
  int    *membership;    /* [numObjs] */
  float **objects;       /* [numObjs][numCoords] data objects */
  float **clusters;      /* [numClusters][numCoords] cluster center */
  float   threshold;
  int     loop_iterations;

  
  /* some default values */
  threshold        = 0.001;
  numClusters      = 0;
  isBinaryFile     = 0;
  is_output_timing = 0;

  // user defined
  char filename[1024] = "./modules/kmeans-master/Image_data/edge100.txt";
  numClusters = 10;

  /* read data points from file ------------------------------------------*/
  objects = file_read(isBinaryFile, filename, &numObjs, &numCoords);
  if (objects == NULL) exit(1);

  /* start the timer for the core computation -----------------------------*/
  /* membership: the cluster id for each data object */
  membership = (int*) malloc(numObjs * sizeof(int));
  assert(membership != NULL);

  clusters = seq_kmeans(objects, numCoords, numObjs, numClusters, threshold,
                        membership, &loop_iterations);

  free(objects[0]);
  free(objects);

  /* output: the coordinates of the cluster centres ----------------------*/
  file_write(filename, numClusters, numObjs, numCoords, clusters,
             membership);

  free(membership);
  free(clusters[0]);
  free(clusters);
}
