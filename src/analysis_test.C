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

using namespace std;

TEST(Stats, pearson)
{
  float array[100];
  float array2[100];
  float array3[100];

  for (int i = 0; i < 100; i++) {
    array[i] = (float) i;
    array2[i] = 2. * (float) i;
    array3[i] = (float) i * (float) i;
  }

  ASSERT_TRUE((mean(array, 100) - 49.5) < 0.001);
  ASSERT_TRUE((pearsonr(array, array2, 100) - 1) < 0.01);

  float p = pearsonr(array, array3, 100);
  ASSERT_TRUE((p - 0.967) < 0.01);
}
