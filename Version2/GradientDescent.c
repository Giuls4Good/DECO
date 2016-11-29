//Program the entire LASSO gradient descent routine in C
#include<stdlib.h>    //needed for malloc
#include<stdio.h>     //needed for printf

// @Input: give pointers to the X-matrix, the Y-matrix, and beta_m of beta=(beta_1, beta_2, ... beta_M)
//         where beta_m is the parameter vector corresponding to the m-th partition.
// @Output: Pointer to beta_m, but only to check if this operation makes sense
double * GradientDescentC(double *restrict X, double *restrict Y, double *restrict beta_m,
                          unsigned int n, unsigned int p, unsigned int m, double precision){


  //return the modified beta_m
  return beta_m;
}



/*
int main(void){
  //define some array dimensions and the chunk-number
  unsigned int n=100, p=50, m=1;
  //create arrays to test your functions in the heap
  double *X, *Y, *beta_m, precision = 0.0001;
  Y = malloc(n*sizeof(double));
  X = malloc(n*p*sizeof(double));
  beta_m = malloc(p*sizeof(double));

  //call the function
  GradientDescentC(X, Y, beta_m, n, p, m, precision);

  //print to the console
  printf("\n %f \n", &beta_m[0]);

  //free the memory allocated to those arrays
  free(Y); free(X); free(beta_m);

  return 0;
}
*/
