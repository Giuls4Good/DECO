//Program the entire LASSO gradient descent routine in C
#include<stdlib.h>    //needed for malloc
#include<stdio.h>     //needed for printf
#include<stdbool.h>   //needed for boolean values
#include<cblas.h>     //needed for CBLAS
//#include<lapacke.h>   //needed for LAPACKE (interface with fortran)
//#include<stddef.h>    //needed for standard definitions

//global variable
const int GD_threshold = 200;    //threshold for analytical vs Gradient Descent OLS
const int max_iter = 1000;       //maximum iterations in the inner loop of CoordinateDescentC
//const double alpha = 1;             //alpha for the function calls of dgemm, dgemv
//const double beta = 0;              //beta for the function calls of dgemm, dgemv


//just some test function
double testFun(const double *restrict X){
  printf("\ntestFun: This function has been called\n");
  return 0.1;
}

//just a printer function for my matrices/vectors
void myPrinter(double * object, int dim1, int dim2){
  printf("\n");
  for(int i = 0; i<dim1; ++i){
    for(int j=0; j<dim2; ++j){
      printf("%f ", *(object+i+j));
    }
    printf("\n");
  }
}

// @Input: give pointers to the X-matrix, the Y-matrix, and beta_m of beta=(beta_1, beta_2, ... beta_M)
//         where beta_m is the parameter vector corresponding to the m-th partition.
//         the X-columns and Y should be standardized, i.e. should have mean 0 and variance 1
// @Output: Pointer to beta_m, but only to check if this operation makes sense
double * CoordinateDescentC(double *restrict X, double *restrict Y, double *restrict beta_m,
                          const int n, const int p, const int m, const double precision){

  //STEP 0: PREPARATIONS (you need some helpers to perform the operations. These are declared here, because
  //                      you will need a separate one for each chunk of the regressor matrix that is being processed)
  double *OLSHelperMatrix, *OLSHelperVector;
  OLSHelperMatrix = malloc(p*p*sizeof(double));
  OLSHelperVector = malloc(p*sizeof(double));

  //STEP 1: INITIALIZE (use OLS is p<n with Gradient Descent if n>GD_threshold, and with solve if n<thresh)

  if(p<n && n>GD_threshold){          //if p<n and the data is large enough to get computational gains from Gradient
                                      //descent for the OLS solution

    //GRADIENT DESCENT OLS
    //GradientDescentC function here

  }else if(p<n && n<GD_threshold){    //if p<n and the data is small enough to use analytically exact OLS


    //ANALYTICAL OLS
    //Perform BLAS operations to obtain OLS solution here.
    //Three steps: (1) Store X'X (dgemm), (2) Store X'Y (dgemv), (3) Store (X'X)^-1(X'Y) (dgesvx)

   /*
    SYNTAX:   C <- alpha*op(A)%*%op(B) + beta*C; op() can be transpose op(x) = x', hermitian op(x) = x^H, or op(x) = x.
              -CBLAS_LAYOUT: row- or column major? (I.e., row- or columnwise reading?)
              -CBLAS_TRANSPOSE: do we want to multiply transpose of A (B) or not?
              -M, N, K: op(A) is MxK, op(B) is KxN, and C is MxN
              -alpha: what we multiply op(A)%*%op(B) with
              -beta: what we multiply C with
              -(A, lda), (B, ldb), (C, ldc): matrices and their respective row (or column) dimension. Row dimension
                if we have column-major, and column dimension if we have row-major

    http://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3_gaeda3cbd99c8fb834a60a6412878226e1.html#gaeda3cbd99c8fb834a60a6412878226e1

    cblas_dgemm (const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB,
                 const int M, const int N, const int K,
                 const double alpha, const double *A, const int lda, const double *B, const int ldb,
                 const double beta, double *OLSHelperMatrix, const int ldc)

    SYNTAX:   Y := alpha*A%*%X + beta*Y OR Y := alpha*A'%*%x + beta*Y
              -incX: specifies the increment in the elements of X (i.e., if we want to move forward one, two,
                or any number of steps between each computation. E.g., if you only want even-numbered elements
                and the 0-th entry, you could choose incX=2. If you want all elements, incX=1)
              -M, N: A is MxN

    http://www.netlib.org/lapack/explore-html/d7/d15/group__double__blas__level2_gadd421a107a488d524859b4a64c1901a9.html#gadd421a107a488d524859b4a64c1901a9

    cblas_dgemv( 	CBLAS_LAYOUT  	layout, CBLAS_TRANSPOSE  	TransA,
                  const int  	M, const int  	N,
                  const double  	alpha, const double *  	A, const int  	lda,
                  const double *  	X, const int  	incX, const double  	beta,
                  double *  	Y, const int  	incY)

    SYNTAX: Compute solution of Linear System A%*%X = B (solution in X)
            -matrix_layout: is either LAPACK_COL_MAJOR or LAPACK_ROW_MAJOR
            -fact: determines if we want to equilibrate (fact = 'E') or not (fact = 'N').
              => if fact='E', the matrices A and B are overwritten with equilibrated versions (!) depending
                 on the value taken by trans
              => if fact='N', they are not overwritten and nor equilibrated
              => in either case,  the LU decomposition is used to split A into  A = P%*%L%*%U
            -trans: determines if you want a transpose, a hermitian, or just the regular matrix for op(A).
              If trans='N', you get op(A) = A. If trans='T', op(A)=A', and if trans='H', op(A)=A^H
            -n: gives order of matrix A/number of linear equations to solve
            -nrhs: gives the 'number of right hand sides' (nrhs), i.e., the column number of B and X (for OLS: nrhs=1)
            -a: the matrix A, with dimension (lda, n). If equed = 'N', or if fact='N', then A will be unchanged on exit.
              If equed = 'N' and fact = 'E', A will also be unchanged. Otherwise, it will be changed to equilibrated version
            -lda: leading dimension of a (A)
            -af: returns factorization arguments from LU-factorization depending on the values of fact and equed
              In any case, it has to have dimension (ldaf, n)
            -ldaf: leading dimension of af, ldaf >= max(1,n)
            -ipiv: integer array that contains output depending on fact (not interesting here)
            -equed: specified which equlibration was performed.
              => equed='N'    no equilibration
              => equed='R'    row equilibration
              => equed='C'    col equilibration
              => equed='B'    both equilibrated
              Note that equed is only an input argument if fact='F', and otherwise it is an output argument
            -r: output argument if fact != 'F', not important here, dimension n
            -c: output argument if fact != 'F', not important here, dimension n
            -b: gives the Right hand side matrix/vector B. Will be overwritten again unless
              equed='N', dimension (ldb, nrhs)
            -ldb: leading dimension of B again
            -x: dimension (ldx, nrhs), returns the solution to the system
            -ldx:leading dimension of x
            -rcond: output parameter giving reciprocal numbers of the matrix (not needed)
            -ferr: forwards error solution bound on each solution (dimension nrhs) (not needed)
            -berr: backwards error solution bound on each solution (dimension nrhs) (not needed)
            -rpivot: another output parameter that is not needed, dimension 4*n

    http://www.netlib.org/lapack/explore-html/d1/d7e/group__double_g_esing_ga325ea9df9f799b0111549999d01cdb53.html#ga325ea9df9f799b0111549999d01cdb53
    http://www.netlib.org/lapack/explore-html/d4/d48/lapacke__dgesvx_8c_source.html
    RETURN: 0 if operation successful, -i if the i-th argument had an illegal value, i and i<=N: U(i,i)
            is exactly zero. i and i = N+1: U nonsingular, but rcond less than machine precision (so the matrix
            is singular to working precision)

    lapack_int LAPACKE_dgesvx( int matrix_layout, char fact, char trans,
                              lapack_int n, lapack_int nrhs,
                              double* a, lapack_int lda,
                              double* af, lapack_int ldaf,
                              lapack_int* ipiv, char* equed, double* r, double* c,
                              double* b, lapack_int ldb,
                              double* x, lapack_int ldx,
                              double* rcond, double* ferr, double* berr, double* rpivot )
    */

    /*
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, //defined in cblas.h, shortcuts to an integer number
                 p, p, n,                               //notice that we give dimensions for AFTER op() is applied!
                 1, X, n, X, n,                        //unclear if lda, ldb are transpose-dependent. (If so, lda = p instead of n here)
                 0, OLSHelperMatrix, p);
    cblas_dgemv( 	CblasColMajor, CblasTrans,
                  n, p,
                  1, X, n,
                  Y, 1, 0,
                  OLSHelperVector, 1);
    //declare some more variables that will be output but that we do not actually care about.
    //notice: might be problematic that they do not have the right dimension?

    double *af, *r, *c, *rcond, *ferr, *berr, *rpivot;
    lapack_int *ipiv, myInfo, ldaf;

    myInfo = LAPACKE_dgesvx(  LAPACK_COL_MAJOR, 'N', 'N',
                               (lapack_int) p, (lapack_int) 1,
                               OLSHelperMatrix, (lapack_int) p,
                               af, lapack_int ldaf,
                               lapack_int* ipiv, 'N', r, c,
                               OLSHelperVector, (lapack_int) p,
                               beta_m, (lapack_int) p,
                              rcond,  ferr,  berr,  rpivot );
    */

    //print the matrix to check if it works

  }else{                              //if p>n, you have to initialize some other way (randomly, e.g.)//check if your convergence criterion is satisfied

    //randomly initialize here

  }

  bool converged = false;             //initialize the convergence criterion to be not achieved yet
  unsigned int count = 0;             //get a counter for the following while loop

  //STEP 2: LOOP
  while(!converged && count<=max_iter){

    //perform UpdateC operation here

    //perform convergence criterion check here

    count++;                          //increment the loop count by 1
  }


  //free memory of OLSHelpers
  free(OLSHelperMatrix); free(OLSHelperVector);

  //return the modified beta_m
  return beta_m;
}




/*int main(void){
  //define some array dimensions and the chunk-number
  const int n=100, p=50, m=1;
  //create arrays to test your functions in the heap
  double *X, *Y, *beta_m;
  const double precision = 0.0001;
  Y = malloc(n*sizeof(double));
  X = malloc(n*p*sizeof(double));
  beta_m = malloc(p*sizeof(double));

  //check printer function
  myPrinter(beta_m, p, 1);

  //call the function
  CoordinateDescentC(X, Y, beta_m, n, p, m, precision);

  //print to the console
  printf("\n %f \n", &beta_m[0]);
  printf("\n GD threshold is %i \n", GD_threshold);

  //free the memory allocated to those arrays
  free(Y); free(X); free(beta_m);

  return 0;
}*/
