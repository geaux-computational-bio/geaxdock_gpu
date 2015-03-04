#include <assert.h>
#include <math.h>
#include "stats.h"


float sum(float * x, int n)
{
  float s = 0.;
  for (int i = 0; i < n; i++)
    s += x[i];
  return s;
}


float mean(float * x, int n)
{
  float s = sum(x, n);
  assert (n > 0);
  return s / (float) n;
}

float pearsonr(float * x, float * y, int n)
{
  float mean_x = mean(x, n);
  float mean_y = mean(y, n);

  float a = 0., b = 0., c = 0.;
  for (int i = 0; i < n; i++) {
    a += (x[i] - mean_x) * (x[i] - mean_x);
    b += (y[i] - mean_y) * (y[i] - mean_y);
    c += (x[i] - mean_x) * (y[i] - mean_y);
  }

  float p_val = c / (sqrt(a) * sqrt(b));

  return p_val;
}
