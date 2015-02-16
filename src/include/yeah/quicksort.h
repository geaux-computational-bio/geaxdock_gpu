//#ifndef _QUICKSORT_H_
//#define _QUICKSORT_H_

#include <stdio.h>
#include <stdlib.h>

void
Swap (int *x, int *y)
{
  int temp;
  temp = *x;
  *x = *y;
  *y = temp;
}



void
QuickSort (int *array, int *order, int first, int last)
{
  if (first < last) {
    int l = first + 1;
    int r = last;
    int pivot = array[first];

    while (l < r) {
      if (array[l] <= pivot)
        l++;
      else if (array[r] > pivot)
        r--;
      else {
        Swap (&array[l], &array[r]);
        Swap (&order[l], &order[r]);
      }
    }
    if (array[l] < pivot) {
      Swap (&array[l], &array[first]);
      Swap (&order[l], &order[first]);
      l--;
    }
    else {
      l--;
      Swap (&array[l], &array[first]);
      Swap (&order[l], &order[first]);
    }

    QuickSort (array, order, first, l);
    QuickSort (array, order, r, last);
  }
}


//#endif
