#include <stdio.h>
#include <stdlib.h>

__declspec(noinline) void vector_add(float* A, float* B, float* C, int n)
{
   int i;
    __assume_aligned(A,64);
    __assume_aligned(B,64);
    __assume_aligned(C,64);
   for (i=0; i<n; ++i){
      C[i] = A[i] + B[i];
   }
}

int main() {
  const int n=512;
  int i;
  float* A = (float*) _mm_malloc(n*sizeof(float),64);
  float* B = (float*) _mm_malloc(n*sizeof(float),64);
  float* C = (float*) _mm_malloc(n*sizeof(float),64);

  // Initialise the data structures
  for (i=0; i<n; i++){
    A[i] = (float) i;
    B[i] = A[i];
  }

  vector_add(A, B, C, n);

  // Output
  for (i=0; i<n; i++)
    printf("%f %f %f\n", A[i], B[i], C[i]);

  _mm_free(A);
  _mm_free(B);
  _mm_free(C);
}
