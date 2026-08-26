#include <stdio.h>
int main()
{
  int proc=0;
#pragma omp parallel
  proc++;
  printf("ParallelThreads: %i\n",proc);
}
