#include <stdio.h>
#include <string>
#include <vector>
#include <cstring>
#include <cmath>

#include "size.h"
#include "dock.h"
#include "load.h"
#include "util.h"
#include "hdf5io.h"
#include "hdf5io.h"
#include "stats.h"
#include "gtest/gtest.h"
#include "gtest/internal/gtest-internal.h"

extern "C" {
#include "./modules/cluster-1.52a/src/cluster.h" /* The C Clustering Library */
}

using namespace std;
/* ========================================================================= */

void show_data(int nrows, int ncols, double** data, int** mask)
/* Print the data matrix */
{ int i, j;
  printf("============== The gene expression data matrix ================\n");
  for (j = 0; j < ncols; j++) printf("\tCol %d", j);
  printf ("\n");
  for (i = 0; i < nrows; i++)
    { printf("Row %d", i);
      for (j = 0; j < ncols; j++)
        { if (mask[i][j]) printf("\t%5.2g",data[i][j]);
          else printf("\t"); /* mask[i][j]==0, so this data point is missing */
        }
      printf("\n");
    }
  printf("\n");
  return;
}

void
example_hierarchical(int nrows, int ncols, double** data, int** mask,
                     double** distmatrix)
/* Perform hierarchical clustering on genes */
{
  int i;
  const int nnodes = nrows -1;
  double* weight = (double*) malloc(ncols*sizeof(double));
  int* clusterid;
  Node* tree;
  for (i = 0; i < ncols; i++) weight[i] = 1.0;
  printf("\n");
  printf("================ Pairwise average linkage clustering ============\n");
  tree = treecluster(nrows, ncols, data, mask, weight, 0, 'e', 'a', 0);
  if (!tree)
    {
      printf ("treecluster routine failed due to insufficient memory\n");
      free(weight);
      return;
    }
  printf("Node     Item 1   Item 2    Distance\n");
  for ( i = 0; i < nnodes; i++)
    printf("%3d:%9d%9d      %g\n",
           -i-1, tree[i].left, tree[i].right, tree[i].distance);
  printf("\n");
  free(tree);
}

TEST(Linkage, cluster)
{
  int i;
  const int nrows = 13;
  const int ncols = 4;
  double** data = (double**) malloc(nrows*sizeof(double*) );
  int** mask = (int**) malloc(nrows*sizeof(int*));
  double** distmatrix;

  for (i = 0; i < nrows; i++)
    { data[i] = (double*) malloc(ncols*sizeof(double));
      mask[i] = (int*) malloc(ncols*sizeof(int));
    }

  /* Test data, roughly distributed in 0-5, 6-8, 9-12 */
  data[ 0][ 0]=0.1; data[ 0][ 1]=0.0; data[ 0][ 2]=9.6; data[ 0][ 3] = 5.6;
  data[ 1][ 0]=1.4; data[ 1][ 1]=1.3; data[ 1][ 2]=0.0; data[ 1][ 3] = 3.8;
  data[ 2][ 0]=1.2; data[ 2][ 1]=2.5; data[ 2][ 2]=0.0; data[ 2][ 3] = 4.8;
  data[ 3][ 0]=2.3; data[ 3][ 1]=1.5; data[ 3][ 2]=9.2; data[ 3][ 3] = 4.3;
  data[ 4][ 0]=1.7; data[ 4][ 1]=0.7; data[ 4][ 2]=9.6; data[ 4][ 3] = 3.4;
  data[ 5][ 0]=0.0; data[ 5][ 1]=3.9; data[ 5][ 2]=9.8; data[ 5][ 3] = 5.1;
  data[ 6][ 0]=6.7; data[ 6][ 1]=3.9; data[ 6][ 2]=5.5; data[ 6][ 3] = 4.8;
  data[ 7][ 0]=0.0; data[ 7][ 1]=6.3; data[ 7][ 2]=5.7; data[ 7][ 3] = 4.3;
  data[ 8][ 0]=5.7; data[ 8][ 1]=6.9; data[ 8][ 2]=5.6; data[ 8][ 3] = 4.3;
  data[ 9][ 0]=0.0; data[ 9][ 1]=2.2; data[ 9][ 2]=5.4; data[ 9][ 3] = 0.0;
  data[10][ 0]=3.8; data[10][ 1]=3.5; data[10][ 2]=5.5; data[10][ 3] = 9.6;
  data[11][ 0]=0.0; data[11][ 1]=2.3; data[11][ 2]=3.6; data[11][ 3] = 8.5;
  data[12][ 0]=4.1; data[12][ 1]=4.5; data[12][ 2]=5.8; data[12][ 3] = 7.6;

  // /* Some data are actually missing */
  mask[ 0][ 0]=1; mask[ 0][ 1]=1; mask[ 0][ 2]=1; mask[ 0][ 3] = 1;
  mask[ 1][ 0]=1; mask[ 1][ 1]=1; mask[ 1][ 2]=0; mask[ 1][ 3] = 1;
  mask[ 2][ 0]=1; mask[ 2][ 1]=1; mask[ 2][ 2]=0; mask[ 2][ 3] = 1;
  mask[ 3][ 0]=1; mask[ 3][ 1]=1; mask[ 3][ 2]=1; mask[ 3][ 3] = 1;
  mask[ 4][ 0]=1; mask[ 4][ 1]=1; mask[ 4][ 2]=1; mask[ 4][ 3] = 1;
  mask[ 5][ 0]=0; mask[ 5][ 1]=1; mask[ 5][ 2]=1; mask[ 5][ 3] = 1;
  mask[ 6][ 0]=1; mask[ 6][ 1]=1; mask[ 6][ 2]=1; mask[ 6][ 3] = 1;
  mask[ 7][ 0]=0; mask[ 7][ 1]=1; mask[ 7][ 2]=1; mask[ 7][ 3] = 1;
  mask[ 8][ 0]=1; mask[ 8][ 1]=1; mask[ 8][ 2]=1; mask[ 8][ 3] = 1;
  mask[ 9][ 0]=1; mask[ 9][ 1]=1; mask[ 9][ 2]=1; mask[ 9][ 3] = 0;
  mask[10][ 0]=1; mask[10][ 1]=1; mask[10][ 2]=1; mask[10][ 3] = 1;
  mask[11][ 0]=0; mask[11][ 1]=1; mask[11][ 2]=1; mask[11][ 3] = 1;
  mask[12][ 0]=1; mask[12][ 1]=1; mask[12][ 2]=1; mask[12][ 3] = 1;

  show_data(nrows, ncols, data, mask);
  example_hierarchical(nrows, ncols, data, mask, distmatrix);
}
