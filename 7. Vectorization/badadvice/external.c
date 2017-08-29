void simple_assign(int n, int* a, int* b){
  int i;
#pragma ivdep
  for (i=0; i<n; i++)
      a[i] = b[i];
}
