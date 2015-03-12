#ifndef __STATS_H_
#define __STATS_H_

#include <vector>

using namespace std;

float sum(float * x, int n);

float mean(float * x, int n);

float pearsonr(float * x, float * y, int n);

float pearsonr(vector < float > x, vector < float > y);

#endif
