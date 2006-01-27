#include <R.h>
#include <Rinternals.h>


void pairwise_sum(double *a, int *na, double *aa){
  int i, j;
  for(i = 0; i < *na; ++i)
    for(j = 0; j < *na; ++j)
      aa[i * *na + j] = a[i] + a[j];
}
