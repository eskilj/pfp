#include <stdio.h>
#include <stdlib.h>

__declspec(noinline) void vector_add(float* A, float* B, float* C, int n)
{
   int i;
   for (i=0; i<n; ++i){
      C[i] = A[i] + B[i];
   }
}

int main() {
  const int n=512;
  int i;
  float* A = (float*) malloc(n*sizeof(float));
  float* B = (float*) malloc(n*sizeof(float));
  float* C = (float*) malloc(n*sizeof(float));

  // Initialise the data structures
  for (i=0; i<n; i++){
    A[i] = (float) i;
    B[i] = A[i];
  }

  vector_add(A, B, C, n);

  // Output
  for (i=0; i<n; i++)
    printf("%f %f %f\n", A[i], B[i], C[i]);

  free(A);
  free(B);
  free(C);
}