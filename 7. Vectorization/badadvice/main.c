#include <stdio.h>

void my_simple_add(int, int*, int*);

int main(){
  const int n=16;
  int i;
  int A[n];
  int B[n];

  // Initialise the data structures
  for (i=0; i<n; i++){
    B[i]=i;
    A[i]=i;
  }

  simple_assign(n-1, B+1, B);

  // Print the final data
  for (i=0; i<n; i++)
    printf("%2d %2d %2d\n", i, A[i], B[i]);
}
